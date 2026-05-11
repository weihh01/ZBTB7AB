library(tidyverse)
pdf("CTCFmotif_ZBTB7A_ZBTB7B-7A7BCTCF.pdf", width = 2.8, height = 2.8)
# pdf("enrichment_HiC_scale_500bp.pdf", width = 3.2, height = 3.9, bg = "transparent")

# df <- read.table("ChIP_SMC1_C2D5_all_CTCFmotif_Mm_profile.csv", sep = ",", header = F)
# df <- as.data.frame(t(data.frame(df, row.names = 1)))

df <- read.table("CTCFmotif_ZBTB7A_ZBTB7B-7A7BCTCF.tab", sep = "\t", header = F, skip = 2)
df <- as.data.frame(t(data.frame(df[,-2], row.names = 1)))
bins <- seq(from = -1000, to = 990, by = 10)
df <- cbind(bins, df)
colnames(df) <- c("bins", "DMSO","dTAG")

p <- ggplot(data = df, mapping = aes(x = bins)) + 
  geom_line(aes(y = df[,2]/1, color = "DMSO"), size = 0.6) +
  geom_line(aes(y = df[,3]/1, color = "dTAG"), size = 0.6) +
  # geom_line(aes(y = df[,6]/1, color = "Hi-TrAC"), size = 1.2) +
  # geom_line(aes(y = df[,2]/1, color = "WT"), size = 1.2) +
  scale_color_manual(values = c(
    "DMSO" = "darkgoldenrod1",
    "dTAG" = "steelblue"
  )
  )+ 
  guides(color=guide_legend(
                            # labels = c("DNase-C", expression(paste(italic('in situ'), ' Hi-C')), "Micro-C", "BL-Hi-C", "Hi-TrAC"),
                            nrow=5, 
                            byrow=TRUE,
                            # label.hjust = 1,
                            keywidth = 0.45
                            )) +
  xlab("CTCF motif") + 
  ylab("CTCF HiChIP signal") +
  labs(title = "") +
  scale_x_continuous(
    limits = c(-1000, 1000),
    breaks = c(-1000, 0, 1000),
    labels = c("-1kb", "0", "1kb"),
    expand = expansion(mult = c(0, 0.03))
  )+
  scale_y_continuous(
    limits = c(0, 10),
    # sbreaks = seq(0,1,0.2),
    # name = "DNase-C enrichment",
    # sec.axis = sec_axis(~. *5, name="DNase-seq enrichment"),
    expand = c(0,0)
  ) + 
  theme(
    axis.text = element_text(color = "black",size = 14), #face = "bold",
    axis.title = element_text(color = "black",size = 14),
    panel.border = element_blank(), 
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    panel.background = element_blank(), 
    legend.position = c(0.8,0.55),
    rect = element_rect(fill = NA, size=0, color = NA),
    plot.background = element_rect(fill = NA, size=0, color = NA),
    legend.background = element_rect(fill = NA, size=0, color = NA),
    legend.key = element_rect(fill = NA, size=0, color = NA),
    # plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  ) 



p

dev.off()
