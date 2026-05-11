library(VennDiagram)
library(RColorBrewer)
library(eulerr)


pdf('gene_venn.pdf', width=5, height=5)
draw.pairwise.venn(266+414,1282+414,414, c("ZBTB7A","ZBTB7B"), 
                   scaled =T, fill = c("#386bb0", "#b02426"),
                   lty = "blank",cat.fontfamily=c("sans","sans"),
                   fontfamily=rep("sans",3),
                   cex = 2,margin=0.12,
                   cat.default.pos="outer",
                   cat.cex = 2,alpha = 0.7,
                   cat.col = c("#386bb0", "#b02426"))
dev.off()
