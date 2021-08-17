library("ggpubr")

setwd("/Users/Corinna/Documents/PostDoc/Beinart_Lab/Snail_16S_amplicons/Meta-analysis/oligotyping")

data <- read.table("size_vs_richness.txt", header=T)

shapiro.test(data$Size)
shapiro.test(data$Oligotypes)

cor.test(data$Size, data$Oligotypes, method="spearman")

pdf("Size_vs_Richness.pdf")
ggscatter(data, x = "Size", y = "Oligotypes", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman", xlab = "Size (mm)", ylab = "Number of oligotypes")
dev.off()
