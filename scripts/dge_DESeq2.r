library(yaml)
library(DESeq2)
library(tidyverse)

args <- commandArgs(TRUE)
countFilepath <- args[1]
output_dea <- args[2]
output_deg <- args[3]

# load the config file
yaml.file <- yaml.load_file('configs/config.yaml')
control <- yaml.file$CONTROL  # all groups used as control
treat <- yaml.file$TREAT  # all groups used as treat, should correspond to control
meta.file <- yaml.file$SAMPLES
meta.data <- read.csv(meta.file, header = TRUE, sep = '\t')
group.all <- meta.data$condition

DEA <- function(countFilepath,output_dea,output_deg,control, treat) {
    countfile <- read.table(countFilepath,header = TRUE, row.names = 1,sep='\t')
    treat_count <- countfile %>% select(contains(treat))
    control_count <- countfile %>% select(contains(control))
    count.table <- cbind(treat_count,control_count)
    condition <- relevel(factor(group.all[c(which(group.all == treat),which(group.all == control))]), ref = control)
    ## create the DESeqDataSet
    colData = data.frame(row.names = colnames(count.table), condition)
    dds <- DESeqDataSetFromMatrix(count.table, colData = colData, design = ~condition)
    ## Filtering
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    ## perform DEA
    dds <- DESeq(dds)
    ## export the results
    res.dea <- results(dds)
    res.dea <- res.dea[complete.cases(res.dea), ]  # remove any rows with NA
    ##dea
    dea <- as.data.frame(res.dea)
    dea <- dea[order(dea$padj, -abs(dea$log2FoldChange), decreasing = FALSE), ]  # sort the table: ascending of FDR then descending of absolute valued of logFC
    #dea$change[dea$padj>=0.05| dea$padj == "NA" | abs(dea$log2FoldChange) < 0.585] <- "Not" 
    #dea$change[dea$padj < 0.05 & dea$log2FoldChange >= 0.585] <- "Up"
    #dea$change[dea$padj < 0.05 & dea$log2FoldChange < -0.585] <- "Down"
    #dea$change[is.na(dea$change)] <- 'Not'
    write.table(data.frame(gene_id=rownames(dea),dea), output_dea, row.names = F, quote = FALSE, sep = '\t')
    ##deg
    deg <- subset(dea,abs(log2FoldChange) >= 1 & padj < 0.05)
    deg <- deg[order(deg$padj, -abs(deg$log2FoldChange), decreasing = FALSE), ]  # sort the table: ascending of FDR then descending of absolute valued of logFC
    write.table(data.frame(gene_id=rownames(deg),deg), output_deg, row.names = F, quote = FALSE, sep = '\t')
}

# the main function
DEA(countFilepath,output_dea,output_deg,control, treat)
