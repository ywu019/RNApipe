library(limma)
library(yaml)
library(edgeR)

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

DEA <- function(countFilepath, output_dea, output_deg, group.all, control, treat) {
    countfile <- read.table(countFilepath, header = T, row.names = 1)
    count.max <- as.matrix(countfile)
    ## create the dataset
    group <- factor(group.all, levels = c(control, treat))
    design <- model.matrix(~0+group)
    colnames(design)=levels(group)
    rownames(design)=colnames(count.max)
    ## Filtering
    dge <- DGEList(counts=count.max, group = group)
    dge <- calcNormFactors(dge)
    v <- voom(dge,design, normalize="quantile")
    ## perform DEA
    fit <- lmFit(v, design)
    constrasts = paste(rev(levels(group)),collapse = "-") 
    cont.matrix <- makeContrasts(contrasts=constrasts,levels  = design) 
    fit2=contrasts.fit(fit,cont.matrix)
    fit2=eBayes(fit2)
    ## dea
    dea <- topTable(fit2, coef=constrasts, n=Inf)
    dea <- na.omit(dea)
    dea <- data.frame(gene_id=rownames(dea), dea)
    write.table(dea, output_dea, row.names = F, quote = FALSE, sep = '\t')
    ## deg
    deg <- subset(dea, abs(logFC) >= 0.585 & adj.P.Val < 0.05)  #adj.P.Val P.Value
    write.table(deg, output_deg, row.names = F, quote = FALSE, sep = '\t')
}

# the main function
DEA(countFilepath, output_dea, output_deg, group.all, control, treat)

