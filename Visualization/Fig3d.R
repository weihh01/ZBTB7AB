library(ggplot2)
library(openxlsx)
library(reshape2)
aaa <- read.xlsx("loop.xlsx")

colnames(aaa) <- c("type", "ratio", "dataset")
aaa$type <- factor(aaa$type, levels = c("PP","EE","EP","struc"), labels = c("Pro-Pro","Enh-Enh","Enh-Pro","Struc-Struc"))
aaa$dataset <- factor(aaa$dataset, levels = c("up","none"), labels = c("strengthened loops","satble loops"))

ggplot(data = aaa, mapping = aes(x = dataset, y = ratio, fill = type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.85) +
  scale_fill_manual(name = "", values = c("#F1B656", "#397FC7", "#6F6F6F","#4CAF50")) +
  labs(x = "", y = "% Loop tpyes", title = NULL) +
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
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.9), breaks = seq(0, 0.9, 0.3), labels = seq(0, 0.9, 0.3) * 100) +
  theme(legend.position = "right", legend.title = element_text(size = 14), legend.text = element_text(size = 13))

ggsave(filename = "EP.pdf", width = 4.9, height = 3.5, dpi = 600)
