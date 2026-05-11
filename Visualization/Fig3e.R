library("openxlsx")
library(ggsignif)

library(ggplot2)
library(reshape2)

df = read.table('genecount.bed',header=F)
colnames(df) <- c("value","variable")
pdf('genecount2.pdf', width=5, height=5)


ggplot(df, aes(x = variable, y = value)) +
  geom_boxplot(aes(fill = variable), show.legend = FALSE, na.rm = TRUE,outlier.shape = NA) +
  scale_fill_manual(
    values = c("up" = "grey", "none" = "#A67BB0")
  ) +
  scale_x_discrete(
    labels = c("up" = "up", "none" = "none")
  ) +
  theme_bw()+
  #ggtitle("DamID-seq signal")+
  theme(panel.grid=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(panel.border = element_blank())+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        #panel.border=element_rect(size=1.2),                  #???呖??哟?
        axis.text.x = element_text(angle=0,vjust = 1,hjust =0.5,color = "black"),
        axis.text.y = element_text(angle=90,vjust = 1,hjust =0.5,color = "black"),
        #axis.text.y = element_text(size =10),
        #legend.position = c(0.75,0.9))+
        #legend.position = c(0.5,0.9),legend.direction="horizontal")+
        legend.position = "none")+
  theme(legend.title=element_blank())+
  labs(x=NULL,y="Gene count per bp")+
  scale_x_discrete(limits = c("up", "none"))+
  theme(text = element_text(size = 20))+
  stat_compare_means(comparisons = list(c("up", "none")), label.y = c(0.4), method = "wilcox.test", label = "p.signif", na.rm = TRUE, paired = FALSE, alternative = "two.sided")


dev.off()
