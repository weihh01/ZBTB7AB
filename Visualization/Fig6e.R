library("openxlsx")
library(ggsignif)

library(ggplot2)
library(reshape2)
library(ggpubr)

library("viridis")
library(ggpointdensity)

gene_diff <- read.table("TSS.bed",header=F,fill=T)
colnames(gene_diff) <- c("DMSO","dTAG","threshold")

gene_diff$DMSO <- log(gene_diff$DMSO,2)
gene_diff$dTAG <- log(gene_diff$dTAG,2)
library(ggplot2)
pdf('TSS_DMSO_dTAG.pdf', width=4, height=4)
ggplot(gene_diff,aes(x=DMSO,y=dTAG,color=threshold))+
  geom_point(size=0.3,alpha=0.5)+
  theme_bw() +
  scale_color_manual(values = c("7A"="#903529", 
                                "7A7B"="black", 
                                "7B"="#AEACAB"))+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 20),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = NA),
        legend.position = c(0.85, 0.25)   
  ) +
  labs(x=expression(log[2]("DMSO")),
       y=expression(log[2]("dTAG"))) +
  scale_y_continuous(breaks = seq(-8,8,4), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(-8,8,4), expand = c(0, 0))
dev.off() 
