library(ggplot2)

gene_diff <- read.csv("DMSO_dTAG.csv", header=TRUE, quote="\t")
result <- gene_diff
result$threshold <- factor(ifelse(result$pvalue < 0.05 & abs(result$log2FoldChange) > 0.6, ifelse(result$log2FoldChange < -0.6, 'Down', 'Up'), 'NoSignifi'), levels = c('Down', 'NoSignifi', 'Up'))

ggplot(result, aes(x = log2FoldChange, y = -log(pvalue, 10), color = threshold)) +
  geom_point(size = 0.1, alpha = 1, shape = 16) +
  scale_color_manual(values = c("#000000", "#AEACAB", "#903529")) +
  theme_bw() +
  theme_classic() +
  annotate("text", x = 0.5, y = 4.5, hjust = 1, label = paste0(format(nrow(result[result$threshold == "Down",]), big.mark = ","), " sites"), size = 6.5, color = "black") +
  annotate("text", x = 0.5, y = 4, hjust = 1, label = paste0(format(nrow(result[result$threshold == "NoSignifi",]), big.mark = ","), " sites"), size = 6.5, color = "#AEACAB") +
  annotate("text", x = 0.5, y = 3.5, hjust = 1, label = paste0(format(nrow(result[result$threshold == "Up",]), big.mark = ","), " sites"), size = 6.5, color = "#903529") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.line.x = element_line(linetype = 1, color = "black", size = 0),
        axis.line.y = element_line(linetype = 1, color = "black", size = 0),
        axis.ticks.x = element_line(color = "black", size = 0.5),
        axis.ticks.y = element_line(color = "black", size = 0.5),
        axis.title.x = element_text(size = 20, vjust = 0),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, colour = 'black'),
        axis.text.y = element_text(size = 20, colour = 'black'),
        plot.title = element_text(colour = "black", size = 20, hjust = 0.5),
        plot.margin = unit(c(0.5, 0.3, 0.3, 0.3), "lines"),
        legend.position = "none") +
  labs(title = "Pol II ChIP-seq", x = expression("log"[2] * " (dTAG / DMSO)"), y = expression("-log"[10] * " (pvalue)"), size = 20) +
  scale_x_continuous(breaks = seq(-3, 3, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 6, 1), expand = c(0, 0))+
  coord_fixed(ratio = 6 / 6 * 1.0, xlim = c(-3, 3), ylim = c(0, 6))

ggsave("volcano.pdf", width = 5, height = 5)
