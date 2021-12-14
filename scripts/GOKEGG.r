library(ggplot2)
library(clusterProfiler)
library(tidyverse)

args <- commandArgs(TRUE)
data <- args[1]
#GO_txt <- args[2]
#KEGG_txt <- args[3]
KEGG_plot <- args[2]
GO_plot <- args[3]
GO_barplot <- args[4]

######### seperate GO and KEGG
kobas <- read.delim(data, comment.char = '#',sep = "\t")
colnames(kobas)=c("Term","Database","ID","Input_num","Background_num","P_Value","Corrected_P_Value","Input","Hyperlink")
go <- subset(kobas, kobas$Database == "Gene Ontology")
KEGG <- subset(kobas, kobas$Database == "KEGG PATHWAY")
GOid <- as.matrix(go[,3])
ont <- go2ont(GOid)
GO <- merge(ont,go,by.x = 1,by.y = 3)
#save
#write.table(GO,GO_txt,sep = '\t',quote = F,row.names = F)
#write.table(KEGG,KEGG_txt,sep = '\t',quote = F,row.names = F)
##################################################### KEGG plot
KEGG <- KEGG[order(KEGG$P_Value),]
KEGG <- KEGG[which(KEGG$P_Value<0.05),]
richFactor <- as.numeric(KEGG$Input_num) / as.numeric(KEGG$Background_num)
KEGG_data <- cbind(KEGG, richFactor)
KEGG_data <- KEGG_data[order(KEGG_data$richFactor),]
KEGG_top <- subset(head(KEGG_data,30))
KEGG_draw <- KEGG_top[,c(1,4,6,7,10)]
KEGG_draw$Term <- factor(KEGG_draw$Term,
                         levels = as.character(KEGG_draw$Term))
KEGG.plot <- ggplot(data = KEGG_draw,aes(x = richFactor, y = Term)) +
  geom_point(aes(size = as.numeric(Input_num),
                 color = as.numeric(P_Value)))+
  scale_color_gradient(low="blue",high="red") +
  labs(x="Gene Ratio",
       #y="Pathway name",
       y="",
       color="p.adjust",
       size="Count") +
  theme_bw() +
  theme(axis.text = element_text(size=12, colour = 'black'),
        axis.title = element_text(size = 12, colour = 'black'),
        panel.grid.major.x = element_line(linetype = 1),
        panel.grid.minor.x = element_line(linetype = 1),
        panel.grid.major.y = element_line(linetype = 1),
        panel.grid.minor.y = element_line(linetype = 1))
ggsave(KEGG.plot, width = 9, height = 6,filename = KEGG_plot,device = "pdf")

# draw the GO plot
GO_sort <- GO[order(GO$P_Value),]
GO_sort_p <- GO_sort[which(GO_sort$P_Value<0.05),]
richFactor <- as.numeric(GO_sort_p$Input_num) / as.numeric(GO_sort_p$Background_num)
GO_data <- cbind(GO_sort_p, richFactor)
GO_data <- GO_data[order(GO_data$richFactor),]
GO_top <- subset(head(GO_data,30))
GO_draw <-GO_top[,c(1,3,5,6,7,8,11)]
GO_draw$Term <- factor(GO_draw$Term,
                       levels = as.character(GO_draw$Term))

# open the figure frame
GO.plot <- ggplot(data = GO_draw,
                  aes(x = richFactor, y = Term)) +
  geom_point(aes(size = as.numeric(Input_num),
                 color = as.numeric(P_Value)))+
  scale_colour_gradient(low="blue",high="red") + 
  labs(x="Gene Ratio",
       #y="Pathway name",
       y="",
       color="p.adjust",
       size="Count"
       # size="Gene number", title=paste("Top",dim(kegg_draw)[1], "of enrichment pathway", sep = ' ')
  ) +
  theme_bw() + 
  theme(axis.text = element_text(size=12, colour = 'black'),
        axis.title = element_text(size = 12, colour = 'black'),
        panel.grid.major.x = element_line(linetype = 1),
        panel.grid.minor.x = element_line(linetype = 1),
        panel.grid.major.y = element_line(linetype = 1),
        panel.grid.minor.y = element_line(linetype = 1))

ggsave(GO.plot, width = 9, height = 6,filename =GO_plot,device = "pdf")

##################################################### draw the GO barplot
BP <- subset(GO_sort,GO_sort$Ontology == "BP" & GO_sort$P_Value < 0.05)
CC <- subset(GO_sort,GO_sort$Ontology == "CC" & GO_sort$P_Value < 0.05)
MF <- subset(GO_sort,GO_sort$Ontology == "MF" & GO_sort$P_Value < 0.05)
mer <- rbind(head(BP,10),head(CC,10),head(MF,10),make.row.names=F)
dorder = factor(as.integer(rownames(mer)),labels=mer$Term)
GO.barplot <- ggplot(mer,aes(x=Term, y=-log(as.numeric(P_Value)),fill=Ontology)) + 
  geom_bar(stat="identity", position=position_dodge(0.6),
           width=0.5, aes(x=dorder)) +
  coord_flip() +
  labs(x="",y="-log10(Pvalue)")+
  scale_y_log10(breaks = c(0,1,10,100,1000)) +
  scale_fill_discrete(name="Ontology") +
  theme(panel.background = element_rect(fill = "transparent",colour = T),
        axis.text = element_text(size=12, colour = 'black'),
        axis.title = element_text(size = 12, colour = 'black'),
        panel.grid.major.x = element_line(linetype = 1),
        panel.grid.minor.x = element_line(linetype = 1),
        panel.grid.major.y = element_line(linetype = 1),
        panel.grid.minor.y = element_line(linetype = 1))
ggsave(GO.barplot, width = 9, height = 6,filename = GO_barplot,device = "pdf")