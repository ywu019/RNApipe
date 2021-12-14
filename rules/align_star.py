import yaml 
import os
import pandas as pd
import numpy as np

with open("configs/config.yaml") as yamlfile:
    config=yaml.load(yamlfile)

if config["TRIM"] == "fastp":
    path = "/trim/fastp"
elif config["TRIM"] == "cutadapt":
    path = "/trim/cutadapt"

END=config["END"]
output_path = config["OUTPUTPATH"]
samples=pd.read_csv(config["SAMPLES"],dtype=str,sep='\t').set_index(["condition","replicate"],drop=False)
condition=np.unique(np.array(pd.read_table(config["SAMPLES"], header = 0)['condition']))
replicate=np.unique(np.array(pd.read_table(config["SAMPLES"], header = 0)['replicate']))
pairend=["R1","R2"]


rule all:
    input:
        bam=expand(output_path + "/align/star/{condition}_{replicate}.sort.bam",condition=condition,replicate=replicate)

rule star_index:
    input:
        genome=config["GENOME"],
        gtf=config["ANNOTATION"]
    output:
        star_index=directory(os.path.dirname(config["GENOME"])+"/star/")
    params:
        index=os.path.dirname(config["GENOME"])+"/star/"
    shell:
        """
        STAR --runThreadN 10 --runMode genomeGenerate --genomeDir {params.index} \
            --genomeFastaFiles {input.genome} --sjdbGTFfile {input.gtf} --genomeSAindexNbases 11
        """

if END=="pair":
    rule star:
        input:
            TrimR1_fastq=lambda wildcards:expand(output_path+path+"/{condition}_{replicate}_R1_trimmed.fastq.gz",condition=wildcards.condition,replicate=wildcards.replicate),
            TrimR2_fastq=lambda wildcards:expand(output_path+path+"/{condition}_{replicate}_R2_trimmed.fastq.gz",condition=wildcards.condition,replicate=wildcards.replicate),
            star_index=os.path.dirname(config["GENOME"])+"/star/"
        output:
            bam=output_path + "/align/star/{condition}_{replicate}.sort.bam"
        params:
            prefix=output_path + "/align/star/{condition}_{replicate}."
        shell:
            """
            STAR --outFilterType BySJout --runThreadN 10 --outFilterMismatchNmax 2 --genomeDir {input.star_index}\
                --readFilesCommand gunzip -c --readFilesIn {input.TrimR1_fastq} {input.TrimR2_fastq} --outFileNamePrefix {params.prefix}\
                    --limitBAMsortRAM 45000000000 --outSAMtype BAM SortedByCoordinate\
                        --quantMode  GeneCounts --outFilterMultimapNmax 1\
                            --outFilterMatchNmin 16
            mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam}
            mv {params.prefix}ReadsPerGene.out.tab {params.prefix}count
            """
else:
    rule star:
        input:
            Trim_fastq=lambda wildcards:expand(output_path+path+"/{condition}_{replicate}_trimmed.fastq.gz",condition=wildcards.condition,replicate=wildcards.replicate),
            star_index=os.path.dirname(config["GENOME"])+"/star/"
        output:
            bam=output_path + "/align/star/{condition}_{replicate}.sort.bam"
        params:
            prefix=output_path + "/align/star/{condition}_{replicate}."
        shell:
            """
            STAR --outFilterType BySJout --runThreadN 10 --outFilterMismatchNmax 2 --genomeDir {input.star_index}\
                --readFilesCommand gunzip -c --readFilesIn {input.Trim_fastq} --outFileNamePrefix {params.prefix}\
                    --limitBAMsortRAM 2000000000 --outSAMtype BAM SortedByCoordinate\
                        --quantMode GeneCounts --outFilterMultimapNmax 1\
                            --outFilterMatchNmin 16
            mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam}
            mv {params.prefix}ReadsPerGene.out.tab {params.prefix}count
            """

