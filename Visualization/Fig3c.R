library(pheatmap)

txt <- readLines("dTAG_big.txt")
txt_clean <- gsub("\\[|\\]", "", txt)
txt_split <- strsplit(txt_clean, ",")
mat1 <- do.call(rbind, lapply(txt_split, function(x) as.numeric(trimws(x))))
pdf('dTAG_big.pdf', width=4, height=4)
rng <- range(mat1, na.rm = TRUE)
pheatmap(
  mat1,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none",
  color = colorRampPalette(c("white", "#fcbba1", "#fb6a4a", "#cb181d"))(100),
  breaks = seq(rng[1], rng[2], length.out = 101),
  border_color = NA,
  legend = TRUE
)
dev.off()



txt <- readLines("DMSO_big.txt")
txt_clean <- gsub("\\[|\\]", "", txt)
txt_split <- strsplit(txt_clean, ",")
mat2 <- do.call(rbind, lapply(txt_split, function(x) as.numeric(trimws(x))))
pdf('DMSO_big.pdf', width=4, height=4)
bk <- c(
  seq(0, 2, length.out = 20),
  seq(2.1, 6, length.out = 20),
  seq(6.1, 9, length.out = 20),
  seq(9.1, 13, length.out = 41)
)

pheatmap(
  mat2,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none",
  color = colorRampPalette(c("white", "#fcbba1", "#fb6a4a", "#cb181d"))(100),
  breaks = seq(rng[1], rng[2], length.out = 101),
  border_color = NA,
  legend = TRUE
)
dev.off() 






txt <- readLines("dTAG_none.txt")
txt_clean <- gsub("\\[|\\]", "", txt)
txt_split <- strsplit(txt_clean, ",")
mat1 <- do.call(rbind, lapply(txt_split, function(x) as.numeric(trimws(x))))
pdf('dTAG_none.pdf', width=4, height=4)


rng <- range(mat1, na.rm = TRUE)

pheatmap(
  mat1,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none",
  color = colorRampPalette(c("white", "#fcbba1", "#fb6a4a", "#cb181d"))(100),
  breaks = seq(rng[1], rng[2], length.out = 101),
  border_color = NA,
  legend = TRUE
)
dev.off()



txt <- readLines("DMSO_none.txt")
txt_clean <- gsub("\\[|\\]", "", txt)
txt_split <- strsplit(txt_clean, ",")
mat2 <- do.call(rbind, lapply(txt_split, function(x) as.numeric(trimws(x))))
pdf('DMSO_none.pdf', width=4, height=4)
bk <- c(
  seq(0, 2, length.out = 20),
  seq(2.1, 6, length.out = 20),
  seq(6.1, 9, length.out = 20),
  seq(9.1, 13, length.out = 41)
)
pheatmap(
  mat2,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none",
  color = colorRampPalette(c("white", "#fcbba1", "#fb6a4a", "#cb181d"))(100),
  breaks = seq(rng[1], rng[2], length.out = 101),
  border_color = NA,
  legend = TRUE
)
dev.off() 


