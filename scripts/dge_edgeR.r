library(edgeR)
library(yaml)
library(tidyverse)

args <- commandArgs(TRUE)
countFilepath <- args[1]
output_dea <- args[2]
output_deg <- args[3]

yaml.file <- yaml.load_file('configs/config.yaml')
control <- yaml.file$CONTROL
treat <- yaml.file$TREAT
meta.file <- yaml.file$SAMPLES
meta.data <- read.csv(meta.file, header = TRUE, sep = '\t')
group.all <- meta.data$condition

DEA <- function(countFilepath, output_dea, output_deg, group.all, control, treat) {
    countfile <- read.table(countFilepath, header = T, row.names = 1)
    count.max <- as.matrix(countfile)
    group <- factor(group.all, levels = c(control, treat))
    design <- model.matrix(~group)
    dgelist <- DGEList(counts = count.max, group = group)
    ## Filtering
    keep <- rowSums(cpm(dgelist) > 1) >= 2
    dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
    dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
    dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
    ## perform DEA
    fit <- glmFit(dge, design, robust = TRUE)
    ## export the results
    lrt <- glmLRT(fit)
    topTags(lrt)
    lrt.all <- topTags(lrt, n = nrow(dgelist$counts))
    ## dea
    dea <- data.frame(gene_id=rownames(lrt.all), lrt.all)
    write.table(dea, output_dea, row.names = F, quote = FALSE, sep = '\t')
    deg <- subset(dea,abs(logFC) >= 1 & FDR < 0.05)
    write.table(deg, output_deg, row.names = F, quote = FALSE, sep = '\t')
}

# the main function
DEA(countFilepath, output_dea, output_deg, group.all, control, treat)
















