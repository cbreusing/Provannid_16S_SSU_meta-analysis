library(ggplot2)
library(plyr)
library(dplyr)
library(devtools)
library(phyloseq)
library(ape)
library("biomformat")
library(vegan)
library(RColorBrewer)
library(circlize)
library(rbiom)

setwd("/Users/Corinna/Documents/PostDoc/Beinart_Lab/16S_amplicons/Meta-analysis/oligotyping")

# Import biom and mapping files
biomFile <- import_biom("oligo-table.biom", parseFunction = parse_taxonomy_default)
mapFile <- import_qiime_sample_data("16S_map.txt")
biomMapFile <- merge_phyloseq(biomFile, mapFile)

# Import tree file
treefile <- read_tree("oligo_tree.nwk")

# Import representative sequences and remove non-OTU information
repsetFile <- Biostrings::readDNAStringSet("oligo_sequences.fasta")
names(repsetFile) <- gsub("\\s.+$", "", names(repsetFile))

biomMapTree <- merge_phyloseq(biomMapFile, treefile)

# Create full phyloseq object
phyloseq <- merge_phyloseq(biomMapTree, repsetFile)
colnames(tax_table(phyloseq)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Symbiont")

symbionts = prune_samples(sample_sums(phyloseq) > 1000, phyloseq)
symbionts = filter_taxa(symbionts, function(x) sum(x) > 0, TRUE)
symbiontstrans <- transform_sample_counts(symbionts, function(x) x/sum(x))

otu <- t(otu_table(symbionts))
spec <- specaccum(otu)

pdf("Specaccum_filtered.pdf")
plot(spec, ci.type="poly", col = "black", lwd=2, ci.lty=0, ci.col="grey", ci = 1.96, xlab = "# individuals", ylab = "# oligotypes")
dev.off()

#Complete dataset
col1 <- colorRampPalette(c("lightcyan2", "cyan4")) (6)
col2 <- "darkslategray"
col3 <- colorRampPalette(c("lightblue1", "midnightblue")) (9)
col4 <- colorRampPalette(c("lavender", "slateblue")) (5)
col5 <- colorRampPalette(c("lemonchiffon", "darkorange")) (19)
col6 <- colorRampPalette(c("mistyrose", "darkred")) (20)

col <- c(col1, col2, col3, col4, col5, col6)

pdf("Fractional_abundance_filtered.pdf", width=30, height=14)
plot_bar(symbiontstrans, x= "Sample", fill = "Symbiont") + facet_grid(. ~ Host, scales = "free", space = "free") + geom_bar(aes(color=Symbiont, fill=Symbiont), stat='identity', position='stack') + scale_color_manual(values = col) + scale_fill_manual(values = col) + ylab("Fractional abundance") + theme_bw() + theme(axis.title = element_text(size=15, face="bold")) + theme(axis.text.x = element_blank(), axis.text.y = element_text(size=13)) + theme(legend.text = element_text(size = 13)) + theme(legend.title = element_text(size = 15, face="bold")) + theme(strip.text.x = element_blank()) + theme(legend.key = element_rect(size = 0.4)) + theme(legend.position="bottom")
dev.off()

colors <- c("pink1", "dodgerblue1", "goldenrod1", "red1", "maroon4", "cyan4", "midnightblue")

# PCoA ordination plots
PCoA <- ordinate(symbiontstrans, "PCoA", "unifrac", weighted = TRUE)
pbc1 <- plot_ordination(symbiontstrans, PCoA)
pdf("PCoA_Unifrac_filtered.pdf")
pbc1 + scale_fill_manual(values = colors, name="Host") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Host)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

boucheti = subset_samples(symbionts, Host == "A. boucheti")
boucheti = filter_taxa(boucheti, function(x) sum(x) > 0, TRUE)
kojimai = subset_samples(symbionts, Host == "A. kojimai")
kojimai = filter_taxa(kojimai, function(x) sum(x) > 0, TRUE)
strummeri = subset_samples(symbionts, Host == "A. strummeri")
strummeri = filter_taxa(strummeri, function(x) sum(x) > 0, TRUE)
hessleri = subset_samples(symbionts, Host == "A. hessleri")
hessleri = filter_taxa(hessleri, function(x) sum(x) > 0, TRUE)
nautilei = subset_samples(symbionts, Host == "I. nautilei")
nautilei = filter_taxa(nautilei, function(x) sum(x) > 0, TRUE)
bouchetitrans <- transform_sample_counts(boucheti, function(x) x/sum(x))
kojimaitrans <- transform_sample_counts(kojimai, function(x) x/sum(x))
strummeritrans <- transform_sample_counts(strummeri, function(x) x/sum(x))
hessleritrans <- transform_sample_counts(hessleri, function(x) x/sum(x))
nautileitrans <- transform_sample_counts(nautilei, function(x) x/sum(x))

colorsb <- c("midnightblue", "goldenrod1", "red2", "cornflowerblue", "lightseagreen", "darkorange")

PCoAb <- ordinate(bouchetitrans, "PCoA", "unifrac", weighted = TRUE)

pbcb <- plot_ordination(bouchetitrans, PCoAb)
pdf("PCoA_Unifrac_boucheti_filtered.pdf")
pbcb + scale_fill_manual(values = colorsb, name="Basin") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Basin)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

colorsk <- c("midnightblue", "goldenrod1", "red2", "cornflowerblue", "maroon4", "darkorange")
colorsk2 <- c("midnightblue", "red2", "goldenrod1")
colorsk3 <- c("midnightblue", "darkorange", "cornflowerblue")

distk <- UniFrac(kojimaitrans, weighted = TRUE)
PCoAk <- cmdscale(distk, eig=TRUE)

var1 <- round(PCoAk$eig[1]/sum(PCoAk$eig), digits=3)*100
var2 <- round(PCoAk$eig[2]/sum(PCoAk$eig), digits=3)*100
x_axis <- paste("Axis.1"," ","[",var1,"%]",sep="")
y_axis <- paste("Axis.2"," ","[",var2,"%]",sep="")

pbck <- plot_ordination(kojimaitrans, PCoAk)
pdf("PCoA_Unifrac_kojimai_filtered.pdf")
pbck + scale_fill_manual(values = colorsk, name="Basin") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Basin)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15)) + labs(x = x_axis, y = y_axis)
dev.off()

ELSCk = subset_samples(kojimai, Basin == "ELSC")
ELSCk = filter_taxa(ELSCk, function(x) sum(x) > 0, TRUE)
ELSCktrans <- transform_sample_counts(ELSCk, function(x) x/sum(x))

distELSCk <- UniFrac(ELSCktrans, weighted = TRUE)
PCoAELSCk <- cmdscale(distELSCk, eig=TRUE)

var1 <- round(PCoAELSCk$eig[1]/sum(PCoAELSCk$eig), digits=3)*100
var2 <- round(PCoAELSCk$eig[2]/sum(PCoAELSCk$eig), digits=3)*100
x_axis <- paste("Axis.1"," ","[",var1,"%]",sep="")
y_axis <- paste("Axis.2"," ","[",var2,"%]",sep="")

pbck2 <- plot_ordination(ELSCktrans, PCoAELSCk)
pdf("PCoA_Unifrac_kojimai_ELSC_filtered.pdf")
pbck2 + scale_fill_manual(values = colorsk2, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15)) + labs(x = x_axis, y = y_axis)
dev.off()

MBk = subset_samples(kojimai, Basin == "Manus Basin")
MBk = filter_taxa(MBk, function(x) sum(x) > 0, TRUE)
MBktrans <- transform_sample_counts(MBk, function(x) x/sum(x))

distMBk <- UniFrac(MBktrans, weighted = TRUE)
PCoAMBk <- cmdscale(distMBk, eig=TRUE)

var1 <- round(PCoAMBk$eig[1]/sum(PCoAMBk$eig), digits=3)*100
var2 <- round(PCoAMBk$eig[2]/sum(PCoAMBk$eig), digits=3)*100
x_axis <- paste("Axis.1"," ","[",var1,"%]",sep="")
y_axis <- paste("Axis.2"," ","[",var2,"%]",sep="")

pbck3 <- plot_ordination(MBktrans, PCoAMBk)
pdf("PCoA_Unifrac_kojimai_MB_filtered.pdf")
pbck3 + scale_fill_manual(values = colorsk3, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15)) + labs(x = x_axis, y = y_axis)
dev.off()

colorss <- c("midnightblue", "goldenrod1", "maroon4")
colorss2 <- c("midnightblue", "red2", "goldenrod1")

PCoAs <- ordinate(strummeritrans, "PCoA", "unifrac", weighted = TRUE)

pbcs <- plot_ordination(strummeritrans, PCoAs)
pdf("PCoA_Unifrac_strummeri_filtered.pdf")
pbcs + scale_fill_manual(values = colorss, name="Basin") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Basin)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

ELSCs = subset_samples(strummeri, Basin == "ELSC")
ELSCs = filter_taxa(ELSCs, function(x) sum(x) > 0, TRUE)
ELSCstrans <- transform_sample_counts(ELSCs, function(x) x/sum(x))

PCoAELSCs <- ordinate(ELSCstrans, "PCoA", "unifrac", weighted = TRUE)

pbcs2 <- plot_ordination(ELSCstrans, PCoAELSCs)
pdf("PCoA_Unifrac_strummeri_ELSC_filtered.pdf")
pbcs2 + scale_fill_manual(values = colorss2, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

colorsh <- c("midnightblue", "goldenrod1", "red2", "cornflowerblue", "lightseagreen")

PCoAh <- ordinate(hessleritrans, "PCoA", "unifrac", weighted = TRUE)

pbch <- plot_ordination(hessleritrans, PCoAh)
pdf("PCoA_Unifrac_hessleri_filtered.pdf")
pbch + scale_fill_manual(values = colorsh, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

colorsn <- c("midnightblue", "lightseagreen")

PCoAn <- ordinate(nautileitrans, "PCoA", "unifrac", weighted = TRUE)

pbcn <- plot_ordination(nautileitrans, PCoAn)
pdf("PCoA_Unifrac_nautilei_filtered.pdf")
pbcn + scale_fill_manual(values = colorsn, name="Basin") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Basin)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

colors2 <- c("royalblue1", "midnightblue", "goldenrod1", "red2", "salmon", "seagreen", "cornflowerblue", "maroon4", "lightseagreen", "darkorange")

# Alpha diversity
pdf("Alpha_diversity_host.pdf")
plot_richness(symbiontstrans, x="Host", measures=c("Shannon", "Simpson")) + scale_fill_manual(values = colors2, name="Basin") + theme_bw() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Basin)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(axis.text.x = element_text(face = "italic"))
dev.off()
