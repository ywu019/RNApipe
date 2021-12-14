library(yaml)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(reshape2)
library(ggrepel)
library(dplyr)
library(cowplot)
library(pheatmap)
library(scales)
library("FactoMineR")
library("factoextra")

args <- commandArgs(TRUE)
file.tpm <- args[1]
output_violin <- args[2]
output_cumulative <- args[3]
output_pca <- args[4]
output_correlation_heatmap <- args[5]
output_fviz <- args[6]
output_tpm_heatmap <- args[7]

yaml.file <- yaml.load_file('configs/config.yaml')
control <- yaml.file$CONTROL
treat <- yaml.file$TREAT
sample <- read.table(yaml.file$SAMPLES, header=T, sep="\t", quote="")
treat_num <- length(which(sample$condition == treat))
control_num <- length(which(sample$condition == control))

plot.cumulative <- function(file.tpm, output_violin, output_cumulative, control, treat) {
    file.tpm <- read.table(file.tpm,header=TRUE,row.names=1)
    pick_row <- apply(file.tpm,1,function(x){
    sum(x<1)<dim(file.tpm)[2]-1
    })
    file.tpm <- file.tpm[pick_row,]
    file.tpm.control <- file.tpm %>% select(contains(control))
    file.tpm.treat <- file.tpm %>% select(contains(treat))
    log2_file.tpm.control <- log2(file.tpm.control+1)
    log2_file.tpm.treat <- log2(file.tpm.treat+1)
    log2_file.tpm.control$control_mean <- apply(log2_file.tpm.control,1,mean)
    log2_file.tpm.treat$treat_mean <- apply(log2_file.tpm.treat,1,mean)
    cumdata1 <- merge(log2_file.tpm.control,log2_file.tpm.treat,by='row.names')
    rownames(cumdata1) <- cumdata1$Row.names
    cumdata1 <- cumdata1 %>% select(treat_mean,control_mean)
    KS = ks.test(x=cumdata1$treat_mean,y=cumdata1$control_mean)
    KS$p.value
    cumdata2 <- melt(cumdata1)
    cumdata2$variable <- factor(cumdata2$variable,levels = c("treat_mean","control_mean"))

    remove_outliers <- function(x, na.rm = TRUE) {
    qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm)
    H <- 1.5 * IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    y
    }
    df2 <- cumdata2 %>% group_by(variable) %>% mutate(value = remove_outliers(value))
    df2 <- df2[complete.cases(df2),]
    
    violin <- ggviolin(df2, x="variable", y="value", fill = "variable",size=0.5,add = "boxplot",add.params = list(fill="white",size=0.5))+ 
    stat_compare_means(comparisons = list(c('treat_mean','control_mean')),label = "p.signif")+
    labs(x='',y="Log2(TPM)")+
    scale_x_discrete(labels = c("treat_mean"=treat,"control_mean"=control))+
    scale_fill_manual(values = c("treat_mean"="#BC3C28","control_mean"="#0072B5"))+
    guides(fill=FALSE)+
    theme(axis.text.x=element_text(size=12),
            plot.margin = margin(0.2, 0.2, 0.01,0.2, "cm"),
            axis.title.y = element_text(size=12),
            panel.grid =element_blank(),
            panel.background = element_rect(color = 'black', fill = 'transparent')
            )
    ggsave(violin,file=output_violin,width = 3, height = 3)


    cumulative <- ggplot(cumdata2,aes(x=value,color=variable))+
    stat_ecdf(size=1)+
    theme(axis.text=element_text(size=12),
            axis.title=element_text(size=12),
            panel.grid =element_blank(),
            panel.background = element_rect(color = 'black', fill = 'transparent'),
            legend.box.just = "right",
            legend.justification = c(0.9,0.85),
            legend.position = c(0.9,0.85),
            legend.key = element_rect(fill = 'transparent',size=0.3))+
    guides(color=guide_legend(title=NULL))+
    scale_y_continuous(limits=c(0,1),breaks = c(0,0.5,1))+
    scale_color_manual(values = c("#BC3C28","#0072B5"),labels=c(treat,control))+
    labs(x="Log2(TPM)",y="Cumulative proportion")
    
    cumulative_plot <- ggdraw(cumulative) + draw_label(paste('P = ',round(KS$p.value,4)), .8, .3,size=12)
    ggsave(cumulative_plot,file=output_cumulative,width = 3, height = 3)
}
plot.cumulative(file.tpm, output_violin, output_cumulative, control, treat)

plot.cor<- function(file.tpm, output_correlation_heatmap) {
    file.tpm <- read.table(file.tpm,header=TRUE,row.names=1)
    pick_row <- apply(file.tpm,1,function(x){
    sum(x<1)<dim(file.tpm)[2]-1
    })
    file.tpm <- file.tpm[pick_row,]
    cormat <- round(cor(file.tpm),2)
    get_upper_tri <- function(cormat){
        cormat[lower.tri(cormat)]<- NA
        return(cormat)}
    upper_tri <- get_upper_tri(cormat)
    melted_cormat <- melt(upper_tri,na.rm = TRUE)

    cor <- ggplot(data=melted_cormat,aes(Var2,Var1,fill=value))+
    geom_tile(color="white")+
    scale_fill_gradient2(low="#3874b1",high = "#bb4b3b",mid="#a98292",midpoint=0.85,space='Lab',name="Cor")+
    coord_fixed()+
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))+
    labs(x="",y="")+
    theme(axis.text=element_text(size=12),
          panel.grid =element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'),
          plot.margin = margin(0.5, 0.01, 0.01,0.01, "cm"),
          legend.text=element_text(size=9),
          legend.background = element_rect(fill = 'transparent'),
          legend.key = element_rect(fill = 'transparent'),
          legend.key.height = unit(0.8,"cm"))
    ggsave(cor,file=output_correlation_heatmap,width = 4, height = 3)
}
plot.cor(file.tpm, output_correlation_heatmap)

plot.pca<- function(file.tpm, output_pca,output_fviz,control_num,treat_num,control,treat) {
    file.tpm <- read.table(file.tpm,header=TRUE,row.names=1)
    pick_row <- apply(file.tpm,1,function(x){
    sum(x<1)<dim(file.tpm)[2]-1
    })
    file.tpm <- file.tpm[pick_row,]
    log2_file.tpm <- log2(file.tpm+1)

    res.pca <- PCA(log2_file.tpm, graph = FALSE,scale = F)
    all_pc <- get_eig(res.pca)
    pc1 <- round(all_pc["Dim.1","variance.percent"])
    pc2 <- round(all_pc["Dim.2","variance.percent"])
    #fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 100))
    all_pc_frame <- as.data.frame(all_pc)
    all_pc_frame$group <- row.names(all_pc_frame) 
    pc_lenth <- dim(all_pc)[1]
    
    fviz_screeplot <- ggplot(all_pc_frame) +
    aes(x=group,y=variance.percent) +
    geom_bar(,stat="identity",position = "dodge",fill="#4682b4") +
    geom_text(aes(label=paste0(round(variance.percent),"%")), 
              position=position_dodge(width = 0.5), vjust=-0.5,size=3)  +
    geom_point(data=all_pc_frame,aes(x=group,y=variance.percent)) + 
    geom_line(data=all_pc_frame,aes(x=group,y=variance.percent,group=1, fill="grey"))+
    theme(axis.text=element_text(size=12),
            axis.title = element_text(size=12),
            panel.grid =element_blank(),
            legend.position = c('none'),
            panel.background = element_rect(color = 'black', fill = 'transparent'),
            )+
    ylim(c(0,100))+labs(x="Dimensions")+
    scale_x_discrete(labels=as.character(c(1:pc_lenth)))

    ggsave(fviz_screeplot,file=output_fviz,width = 3, height = 3)
    ## pca
    sample_mrna <- data.frame(sample=colnames(log2_file.tpm),
                          group=c(rep(treat, treat_num), rep(control, control_num)))
    data <- res.pca[["var"]][["coord"]]
    dat <- merge(data,sample_mrna,by.x='row.names',by.y='sample')
    dat$group <- factor(dat$group,levels=c(treat,control))
    
    pca <- ggplot(dat,aes(x=Dim.1,y=Dim.2))+
    geom_point(size=3,aes(fill=group),shape=21)+
    geom_text_repel(aes(x=Dim.1,y=Dim.2, label=Row.names))+
    scale_fill_manual(values = c("#BC3C28","#0072B5"),labels=c(treat,control))+
    labs(x=paste("PC1 (",pc1,"%)"),y=paste("PC2 (",pc2,"%)"))+
    geom_hline(yintercept = 0,size=0.2)+
    guides(fill=FALSE)+
    theme(axis.text=element_text(size=12),
            axis.title=element_text(size=12),
            panel.grid =element_blank(),
            panel.background = element_rect(color = 'black', fill = 'transparent'),
            legend.key = element_rect(fill = 'transparent'))
    ggsave(pca,file=output_pca,width = 3, height = 3)
}
plot.pca(file.tpm,output_pca,output_fviz,control_num,treat_num,control, treat)

pheatmap.tpm<- function(file.tpm, output_tpm_heatmap) {
    file.tpm <- read.table(file.tpm,header=TRUE,row.names=1)
    file.tpm.control <- file.tpm %>% select(contains(control))
    file.tpm.treat <- file.tpm %>% select(contains(treat))
    num.control <- dim(file.tpm.control)[2]
    num.treat <- dim(file.tpm.treat)[2]
    file.tpm <- cbind(file.tpm.treat,file.tpm.control)
    annotation_col = data.frame(Type = factor(c(rep("treat",num.treat),rep("control", num.control))))
    rownames(annotation_col) = colnames(file.tpm)
    ann_colors = list(Type = c(treat="#BC3C28",control="#0072B5"))
  
    pick_row <- apply(file.tpm,1,function(x){
      sum(x<1)<dim(file.tpm)[2]-1
    })
    file.tpm <- file.tpm[pick_row,]
    log2_file.tpm <- log2(file.tpm+1)
    pdf(file = output_tpm_heatmap, width = 4, height = 6)
    pheatmap(log2_file.tpm,
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             show_rownames = F,show_colnames = T,
             cluster_cols = T,
             cluster_rows=T,
             treeheight_row=8,
             treeheight_col=5,
             cellwidth = 28,border=FALSE,
             color =colorRampPalette(c("#0072B5","white","#BC3C28"))(100),
             fontsize_col=8, 
             legend_breaks = c(-1,0,1),
             scale = "row",angle_col=45, 
             annotation_legend = F,
             clustering_distance_rows = 'euclidean', 
             clustering_method = 'single',
             main = "Heatmap of all expressed genes"
    )
    dev.off()
}

pheatmap.tpm(file.tpm, output_tpm_heatmap)
