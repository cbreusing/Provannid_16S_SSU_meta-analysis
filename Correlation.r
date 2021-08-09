library("ggpubr")

setwd("/Users/Corinna/Documents/PostDoc/Beinart_Lab/16S_amplicons/Meta-analysis/oligotyping")

data <- read.table("depth_vs_richness.txt", header=T)

shapiro.test(data$Depth)
shapiro.test(data$Oligotypes)

cor.test(data$Depth, data$Oligotypes, method="spearman")

pdf("Depth_vs_Richness.pdf")
ggscatter(data, x = "Depth", y = "Oligotypes", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman", xlab = "Depth (m)", ylab = "Number of oligotypes")
dev.off()

data2 <- read.table("size_vs_richness.txt", header=T)

shapiro.test(data2$Size)
shapiro.test(data2$Oligotypes)

cor.test(data2$Size, data2$Oligotypes, method="spearman")

pdf("Size_vs_Richness.pdf")
ggscatter(data2, x = "Size", y = "Oligotypes", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman", xlab = "Size (mm)", ylab = "Number of oligotypes")
dev.off()
