import yaml
import os
import time

with open('configs/config.yaml') as yamlfile:
    config = yaml.load(yamlfile, Loader=yaml.FullLoader)
project = config["PROJECT"]
bigwig = config["BIGWIG"][0]
thread=config["THREAD"]
wait_time=config["TIME"]

# Start the workflow
print("Start RNApipe on project: " + project)

## write the running time in a log file
file_log_time = open("logs/log_running_time.txt", "a+")
file_log_time.write("\nProject name: " + project + "\n")
file_log_time.write("Start time: " + time.ctime() + "\n")
def spend_time(start_time, end_time):
    seconds = end_time - start_time
    hours = seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60
    return "%d:%02d:%02d" % (hours, minutes, seconds)

# Trim
if config["TRIM"] == "fastp":
    print("Start Trimming!")
    start_time = time.time()
    os.system("nice -5 snakemake -j " + thread + " -s rules/qc_fastp.py --latency-wait " + wait_time + " 2>&1 | tee logs/log_fastp.txt")
    end_time = time.time()
    file_log_time.write("Time of running trimming:" + spend_time(start_time, end_time) + "\n")
    print("Trimming is done!")
elif config["TRIM"] == "cutadapt":
    print("Start Trimming!")
    start_time = time.time()
    os.system("nice -5 snakemake -j " + thread + " -s rules/qc_cutadapt.py --latency-wait " + wait_time + " 2>&1 | tee logs/log_cutadapt.txt")
    end_time = time.time()
    file_log_time.write("Time of running trimming:" + spend_time(start_time, end_time) + "\n")
    print("Trimming is done!") 

# Mapping
if config["ALIGN"] == "hisat2":
    print("Start Mapping to genome!")
    start_time = time.time()
    os.system("nice -5 snakemake -j " + thread + " -s rules/align_hisat2.py --latency-wait " + wait_time + " 2>&1 | tee logs/log_hisat2.txt")
    end_time = time.time()
    file_log_time.write("Time of running mapping:" + spend_time(start_time, end_time) + "\n")
    print("Mapping is done!")
elif config["ALIGN"] == "STAR":
    print("Start Mapping to genome!")
    start_time = time.time()
    os.system("nice -5 snakemake -j " + thread + " -s rules/align_star.py --latency-wait " + wait_time + " 2>&1 | tee logs/log_star.txt")
    end_time = time.time()
    file_log_time.write("Time of running mapping:" + spend_time(start_time, end_time) + "\n")
    print("Mapping is done!") 

# bamtobigwig
if config["BIGWIG"] == "yes":
    print("Start transform bam file to bigwig file!")
    start_time = time.time()
    os.system("nice -5 snakemake -j " + thread + " -s rules/bamtobigwig.py --latency-wait " + wait_time + " 2>&1 | tee logs/log_deeptools.txt")
    end_time = time.time()
    file_log_time.write("Time of running:" + spend_time(start_time, end_time) + "\n")
    print("Transforming is done!")

# count
if config["QUANTIFY"] == "featureCounts":
    print("Start Quantifying!")
    start_time = time.time()
    os.system("nice -5 snakemake -j " + thread + " -s rules/featurecounts.py --latency-wait " + wait_time + " 2>&1 | tee logs/log_featurecounts.txt")
    end_time = time.time()
    file_log_time.write("Time of running quantifying:" + spend_time(start_time, end_time) + "\n")
    print("Quantifying is done!")
elif config["QUANTIFY"] == "htseq":
    print("Start Quantifying!")
    start_time = time.time()
    os.system("nice -5 snakemake -j " + thread + " -s rules/htseq.py --latency-wait " + wait_time + " 2>&1 | tee logs/log_htseq.txt")
    end_time = time.time()
    file_log_time.write("Time of running quantifying:" + spend_time(start_time, end_time) + "\n")
    print("Quantifying is done!")

#visualize
print("Start doing visualization!")
start_time = time.time()
os.system("nice -5 snakemake -j " + thread + " -s rules/exp_visualize.py --latency-wait " + wait_time + " 2>&1 | tee logs/log_visualize.txt")
end_time = time.time()
file_log_time.write("Time of running visualization:" + spend_time(start_time, end_time) + "\n")
print("visualization is done!")

# multiqc
print("Start Multiqc!")
start_time = time.time()
os.system("nice -5 snakemake -j " + thread + " -s rules/qualimap.py --latency-wait " + wait_time + " 2>&1 | tee logs/log_qualimap.txt")
end_time = time.time()
file_log_time.write("Time of running multiqc:" + spend_time(start_time, end_time) + "\n")
print("Multiqc is done!")

# dea
if config["DEA"] == "DESeq2":
    print("Start doing DEA!")
    start_time = time.time()
    os.system("nice -5 snakemake -j " + thread + " -s rules/dge_DESeq2.py --latency-wait " + wait_time + " 2>&1 | tee logs/log_dea_DESeq2.txt")
    end_time = time.time()
    file_log_time.write("Time of running DEA genome based:" + spend_time(start_time, end_time) + "\n")
    print("DEA is done!")
elif config["DEA"] == "edgeR":
    print("Start doing DEA!")
    start_time = time.time()
    os.system("nice -5 snakemake -j " + thread + " -s rules/dge_edgeR.py --latency-wait " + wait_time + " 2>&1 | tee logs/log_dea_edgeR.txt")
    end_time = time.time()
    file_log_time.write("Time of running DEA genome based:" + spend_time(start_time, end_time) + "\n")
    print("DEA is done!") 
elif config["DEA"] == "limma":
    print("Start doing DEA!")
    start_time = time.time()
    os.system("nice -5 snakemake -j " + thread + " -s rules/dge_limma.py --latency-wait " + wait_time + " 2>&1 | tee logs/log_dea_limma.txt")
    end_time = time.time()
    file_log_time.write("Time of running DEA genome based:" + spend_time(start_time, end_time) + "\n")
    print("DEA is done!")

#GO KEGG
print("Start doing Functional annotation!")
start_time = time.time()
os.system("nice -5 snakemake --use-conda -j " + thread + " -s rules/go_kegg.py --latency-wait " + wait_time + " 2>&1 | tee logs/log_kobas.txt")
end_time = time.time()
file_log_time.write("Time of running Functional annotation:" + spend_time(start_time, end_time) + "\n")
print("Functional annotation is done!")

#GSEA
print("Start doing GSEA!")
start_time = time.time()
os.system("nice -5 snakemake --use-conda -j " + thread + " -s rules/gsea.py --latency-wait " + wait_time + " 2>&1 | tee logs/log_GSEA.txt")
end_time = time.time()
file_log_time.write("Time of running GSEA:" + spend_time(start_time, end_time) + "\n")
print("GSEA is done!")
