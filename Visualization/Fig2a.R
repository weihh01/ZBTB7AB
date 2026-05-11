library(ggplot2)
library(openxlsx)
library(reshape2)
aaa <- read.xlsx("motif_CTCF.xlsx")

colnames(aaa) <- c("type", "ratio", "dataset")
#aaa$type <- factor(aaa$type, levels = c("CTCFmotif", "noCTCFmotif"), labels = c("CTCFmotif", "noCTCFmotif"))
aaa$dataset <- factor(aaa$dataset, levels = c("ZBTB7A_ZBTB7B", "ZBTB7A_noZBTB7B", "ZBTB7B_noZBTB7A"), labels = c("7A_7B", "7A_no7B", "7B_no7A"))

ggplot(data = aaa, mapping = aes(x = dataset, y = ratio, fill = type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.85) +
  scale_fill_manual(name = "", values = c("#4E79A7","#F28E2B")) +
  labs(x = "", y = "Percentage (%)", title = NULL) +
  guides(fill = guide_legend(ncol = 1)) +
  geom_text(aes(label = paste0(round(ratio * 100, 1))), position = position_stack(vjust = 0.5), size = 4.5, colour = "black") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 0.5),
        axis.line.x = element_blank(),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        axis.text.x = element_text(size = 14, colour = "black", angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, color = "black", hjust = 0.5),
        axis.title.y = element_text(size = 14, color = "black", hjust = 0.5),
        plot.title = element_text(colour = "black", face = "bold",
                                  size = 14, vjust = 2, hjust = 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.4), breaks = seq(0, 0.4, 0.1), labels = seq(0, 0.4, 0.1) * 100) +
  theme(legend.position = "right", legend.title = element_text(size = 14), legend.text = element_text(size = 13))

ggsave(filename = "motif_CTCF.pdf", width = 4.5, height = 3.3, dpi = 600)
