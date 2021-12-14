import yaml
import os
import pandas as pd
import numpy as np

with open("configs/config.yaml") as yamlfile:
    config=yaml.load(yamlfile)

if config["ALIGN"] == "hisat2":
    path = "/align/hisat2"
elif config["ALIGN"] == "STAR":
    path = "/align/star"

END=config["END"]
output_path = config["OUTPUTPATH"]
samples=pd.read_csv(config["SAMPLES"],sep='\t')
condition=np.unique(samples['condition'])
replicate=np.unique(samples['replicate'])
pairend=["R1","R2"]

rule all:
    input:
        tpm_tsv=output_path + "/quantify/htseq/genome_tpm.tsv"

rule resort:
    input:
        sortbam=output_path + path + "/{condition}_{replicate}.sort.bam"
    output:
        resortbam=temp(output_path + path + "/{condition}_{replicate}_name.bam")
    shell:
        """
        samtools sort -@ 10 -n {input.sortbam} -o {output.resortbam}
        """

if END=='pair':
    rule htseq:
        input:
            resortbam=output_path + path + "/{condition}_{replicate}_name.bam",
            gtf=config["ANNOTATION"]
        output:
            summary=temp(output_path + "/quantify/htseq/{condition}_{replicate}.txt")
        shell:
            """
            htseq-count -f bam -r name -s reverse -a 10 -t exon -i gene_id -m union --nonunique=none {input.resortbam} {input.gtf}  > {output.summary}
            """
else:
    rule htseq:
        input:
            resortbam=output_path + path + "/{condition}_{replicate}_name.bam",
            gtf=config["ANNOTATION"]
        output:
            summary=output_path + "/quantify/htseq/{condition}_{replicate}.txt"
        shell:
            """
            htseq-count -f bam -r name -s no -a 10 -t exon -i gene_id -m union --nonunique=none {input.resortbam} {input.gtf} > {output.summary}
            """

rule combine_summary: 
    input: 
        summary=expand(output_path + "/quantify/htseq/{condition}_{replicate}.txt",condition=condition,replicate=replicate)
    output:
        genome_summary = output_path + "/quantify/htseq/genome.tsv.summary",
        genome_count = output_path + "/quantify/htseq/genome_count.tsv"
    params:
        input_path=output_path + "/quantify/htseq/"
    shell: 
        """
        Rscript scripts/combine_count.r {params.input_path} {output.genome_summary}
        sed '/^__/d' {output.genome_summary} > {output.genome_count}
        """

rule get_gene_length:
    input:
        gtf = config["ANNOTATION"],
        genome_count = output_path + "/quantify/htseq/genome_count.tsv"
    output:
        gene_length = temp(output_path + "/quantify/htseq/gene_length.tsv")
    shell:
        """
        python scripts/get_gene_length.py -i {input.gtf} -o {output.gene_length}
        """

rule Counts2tpm:
    input:
        genome_count = output_path + "/quantify/htseq/genome_count.tsv",
        gene_length=output_path + "/quantify/htseq/gene_length.tsv"
    output:
        tpm_tsv=output_path + "/quantify/htseq/genome_tpm.tsv"
    shell:
        """
        python scripts/htseqcounts2tpm.py -i {input.genome_count} -l {input.gene_length} -o {output.tpm_tsv}
        """
