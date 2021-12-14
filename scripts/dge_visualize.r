library(yaml)
library(ggplot2)
library(pheatmap)
library(tidyverse)

args <- commandArgs(TRUE)
file.dea <- args[1] 
file.tpm <- args[2]
output_valcano <- args[3]
output_heatmap <- args[4]

yaml.file <- yaml.load_file('configs/config.yaml')
control <- yaml.file$CONTROL
treat <- yaml.file$TREAT

plot.volcano <- function(file.dea, output_valcano, control, treat) {
  dea.table <- read.table(file.dea, header = TRUE, row.names = 1)
  # volcano plot
  dea.table.volcano <- dea.table  # for better volcano plot, 0 FDRs/padj will be changed to a very low value
  dea.table.volcano["group"] <- ifelse(dea.table.volcano$padj>=0.05| dea.table.volcano$padj %in% NA | abs(dea.table.volcano$log2FoldChange) < 1,'Not',
                             ifelse(dea.table.volcano$padj < 0.05 & dea.table.volcano$log2FoldChange >= 1,'Up','Down'))
  up_num <- length(which(dea.table.volcano["group"]=="Up"))
  down_num <- length(which(dea.table.volcano["group"]=="Down"))
  dea.table.volcano$group <- factor(dea.table.volcano$group,levels=c("Up","Down","Not"))
  # change the 0 padj to a low value (100 times smaller than the minumum non-zero value)
  padj <- dea.table.volcano$padj
  padj.min.non.zero <- min(padj[padj>0])
  dea.table.volcano$padj[padj==0] <- padj.min.non.zero/100
  ## define the range of x-axis and y-axis
  log2FC_lim <- max(abs(dea.table.volcano$log2FoldChange))
  padj_lim <- -log10(min(dea.table.volcano$padj)) # NAs already removed from dea.table

  ggplot2.volcano <- ggplot(dea.table.volcano,aes(log2FoldChange,-1*log10(padj),color = group))+
      geom_point(size=0.8)+
      labs(x="log2(FoldChange)",y="-log10(padj)",title=paste(treat,"vs",control))+ 
      scale_color_manual(values =c('Up'="#BC3C28","Down"="#0072B5","Not"="grey"),
                         labels=c(paste("Up",up_num),paste("Down",down_num),paste("Normal")))+
      geom_hline(yintercept=-log10(0.05),linetype=4)+
      geom_vline(xintercept=c(-1,1),linetype=4)+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=12, hjust = 0.5),
            plot.title=element_text(size=12, hjust = 0.5),
            panel.grid =element_blank(),
            panel.background = element_rect(color = 'black', fill = 'transparent'),
            legend.text=element_text(size=10),
            legend.background = element_rect(fill = 'transparent'),
            legend.key = element_rect(fill = 'transparent'),
            legend.key.height = unit(0.8,"cm"))+
      xlim(-log2FC_lim-1,log2FC_lim+1)+
      guides(color=guide_legend(title=NULL))
  ggsave(ggplot2.volcano, width = 4, height = 3,file = output_valcano)
}
plot.volcano(file.dea, output_valcano, control, treat)


plot.heatmap <- function(file.dea,file.tpm, output_heatmap, control, treat) {
  file.tpm <- read.table(file.tpm,header=TRUE,row.names=1)
  dea.table <- read.table(file.dea, header = TRUE, row.names = 1)
  dea.table <- dea.table[order(dea.table$padj, -abs(dea.table$log2FoldChange), decreasing = FALSE), ]
  gene.id.dea <- row.names(dea.table)
  file.tpm.control <- file.tpm %>% select(contains(control))
  file.tpm.treat <- file.tpm %>% select(contains(treat))
  num.control <- dim(file.tpm.control)[2]
  num.treat <- dim(file.tpm.treat)[2]
  file.tpm <- cbind(file.tpm.treat,file.tpm.control)
  index.deg <- which(row.names(file.tpm) %in% gene.id.dea[1:30])
  file.tpm.deg <- file.tpm[index.deg,]
  gene.id.file.tpm <- rownames(file.tpm.deg)
  annotation_col = data.frame(Type = factor(c(rep("treat",num.treat),rep("control", num.control))))
  rownames(annotation_col) = colnames(file.tpm.deg)
  ann_colors = list(Type = c(treat="#BC3C28",control="#0072B5"))
  pdf(file = output_heatmap, width = 4, height = 5)
  pheatmap(as.matrix(file.tpm.deg),
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           labRow = gene.id.file.tpm,
           fontsize=6,scale="row",
           angle_col = "45",
           treeheight_row=8,treeheight_col=8, fontsize_col=8,
           cellwidth = 28, cellheight = 9,cluster_row = T,border=FALSE,
           legend_breaks = c(-1,0,1),
           main = "Heatmap of the top 30 differential expressed genes",
           annotation_legend = F)
  dev.off()
}
plot.heatmap(file.dea,file.tpm, output_heatmap, control, treat)
