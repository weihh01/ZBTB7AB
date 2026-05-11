library(openxlsx)
library(ggsignif)
library(ggplot2)
library(reshape2)
library(ggpubr)

mymap1 <- mymap2 <- mymap3 <- c("grey", "#A67BB0")
df <- read.table('CTCF_dTAG-DMSO.bed', header = F)
colnames(df) <- c("value", "variable")


ggplot(df, aes(x = variable, y = value)) +
  geom_boxplot(aes(fill = variable), show.legend = F, na.rm = TRUE, outlier.shape = NA) +
  theme_bw() +
  ggtitle("ZBTB7A ChIP-seq signals") +
  scale_fill_manual(values = mymap1) +
  scale_color_manual(values = c("black", "black")) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.8, colour = "black", size = 14),
        axis.text.y = element_text(size = 14, face = "plain", colour = "black", angle = 0),
        axis.title.x = element_text(size = 15, face = "plain", colour = "black", vjust = 0, hjust = 0.5),
        axis.title.y = element_text(size = 15, face = "plain", colour = "black", vjust = 1, hjust = 0.5),
        panel.border = element_blank(),
        axis.line.x = element_line(linetype = 1, color = "black", size = 0.5),
        axis.line.y = element_line(linetype = 1, color = "black", size = 0.5),
        axis.ticks.x = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.5),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines"),
        plot.title = element_text(colour = "black", face = "plain", size = 15, vjust = 0.5, hjust = 0.5))+
  labs(title = "ZBTB7A ChIP-seq", x = NULL, y = expression(log[2]("dTAG/DMSO"))) +
  ylim(-0.3, 0.5)+
  theme(text = element_text(size = 14)) +
  stat_compare_means(comparisons = list(c("up", "none")), label.y = c(0.4), method = "wilcox.test", label = "p.signif", na.rm = TRUE, paired = FALSE, alternative = "two.sided")

ggsave(filename = "CTCF_dTAG-DMSO.pdf", width = 3, height = 3, dpi = 600)
