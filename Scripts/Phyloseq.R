library(ggplot2)
library(plyr)
library(dplyr)
library(devtools)
library(phyloseq)
library(ape)
library(biomformat)
library(vegan)
library(RColorBrewer)
library(circlize)
library(rbiom)
library(forcats)
library(microbiome)
library(LDM)

setwd("/Users/Corinna/Documents/PostDoc/Beinart_Lab/Snail_16S_amplicons/Meta-analysis/oligotyping")

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

# Total sum scaling
symbiontstrans <- transform_sample_counts(symbionts, function(x) x/sum(x))
bouchetitrans <- transform_sample_counts(boucheti, function(x) x/sum(x))
kojimaitrans <- transform_sample_counts(kojimai, function(x) x/sum(x))
strummeritrans <- transform_sample_counts(strummeri, function(x) x/sum(x))
hessleritrans <- transform_sample_counts(hessleri, function(x) x/sum(x))
nautileitrans <- transform_sample_counts(nautilei, function(x) x/sum(x))

# Hellinger transformation
hellinger.symbionts <- transform(symbionts, "hellinger")
hellinger.boucheti <- transform(boucheti, "hellinger")
hellinger.kojimai <- transform(kojimai, "hellinger")
hellinger.strummeri <- transform(strummeri, "hellinger")
hellinger.hessleri <- transform(hessleri, "hellinger")
hellinger.nautilei <- transform(nautilei, "hellinger")

# ASV accumulation curves
otub <- t(otu_table(boucheti))
otuk <- t(otu_table(kojimai))
otus <- t(otu_table(strummeri))
otuh <- t(otu_table(hessleri))
otun <- t(otu_table(nautilei))
otu <- t(otu_table(symbionts))

specb <- specaccum(otub)
speck <- specaccum(otuk)
specs <- specaccum(otus)
spech <- specaccum(otuh)
specn <- specaccum(otun)
spec <- specaccum(otu)

pdf("Specaccum_filtered_combined.pdf")
plot(spec, col = "black", lwd=2, ci = 0, xlab = "# individuals", ylab = "# oligotypes")
plot(specb, col = "dodgerblue1", lwd=2, add = TRUE, ci = 0)
plot(speck, col = "red1", lwd=2, add = TRUE, ci = 0)
plot(specs, col = "cyan4", lwd=2, add = TRUE, ci = 0)
plot(spech, col = "darkorange", lwd=2, add = TRUE, ci = 0)
plot(specn, col = "midnightblue", lwd=2, add = TRUE, ci = 0)
legend("bottomright", legend=c("All species", "A. boucheti", "A. kojimai", "A. strummeri", "A. hessleri", "I. nautilei"), col=c("black", "dodgerblue1", "red1", "cyan4", "darkorange", "midnightblue"), lty=1, cex=1, text.font=c(1, 3, 3, 3, 3, 3), box.lty=0)
dev.off()

#Complete dataset
col1 <- colorRampPalette(c("lightcyan2", "cyan4")) (6)
col2 <- "darkslategray"
col3 <- colorRampPalette(c("lightblue1", "midnightblue")) (9)
col4 <- colorRampPalette(c("lavender", "slateblue")) (5)
col5 <- colorRampPalette(c("lemonchiffon", "darkorange")) (19)
col6 <- colorRampPalette(c("mistyrose", "darkred")) (20)

col <- c(col1, col2, col3, col4, col5, col6)

desired_order <- list("MA7", "MA9", "3", "17", "35", "37", "Al0050", "Al0076", "Al0077", "Al0082", "Al0083", "Al0086", "Al0094", "Al0499", "Al0503", "Al0505", "Al0506", "Al0511", "Al0514", "Al0519", "Al0520", "Al0523", "Al0467", "Al0468", "Al0472", "Al0532", "Al0535", "Al0536", "Al0537", "Al0542", "Al0544", "Al0546", "Al0550", "Al0551", "Al0558", "Al0612", "Al0613", "Al0617", "Al0721", "Al0722", "Al0723", "Al0724", "Al0726", "Al0727", "Al0732", "Al0734", "Al0737", "Al0739", "Al0743", "Al0744", "MB-11", "MB-12", "MB-13", "MB-14", "MB-15", "MB-16", "MB-19", "MB-21", "MB-23", "MB-26", "MB-36", "Al0321", "Al0323", "Al0324", "Al0325", "Al0326", "Al0327", "Al0328", "Al0329", "Al0330", "Al0331", "Al0333", "Al0334", "Al0335", "Al0337", "Al0338", "Al0339", "Al0340", "Al0346", "826", "859", "863", "867", "869", "871", "873", "877", "879", "881", "883", "885", "887", "889", "891", "893", "895", "897", "901", "903", "905", "907", "911", "913", "915", "917", "919", "923", "925", "927", "929", "931", "933", "935", "937", "939", "941", "944", "946", "948", "951", "954", "958", "960", "963", "964", "966", "968", "970", "972", "974", "976", "978", "982", "983", "986", "988", "990", "994", "1097", "1099", "1101", "1103", "1105", "1107", "1109", "1111", "1113", "1115", "1117", "1119", "1121", "1123", "1125", "1127", "1129", "1131", "1135", "1137", "1140", "1147", "1149", "1151", "MA11", "MA13", "MA15", "MA17", "MA172", "MA174", "MA178", "MA182", "MA184", "MA186", "MA188", "MA19", "MA190", "MA192", "MA194", "MA196", "MA198", "MA208", "MA21", "MA210", "MA212", "MA22", "MA23", "MA24", "MA26", "MA28", "MA31", "MA32", "MA34", "MA36", "MA38", "MA40", "MA42", "MA44", "MA46", "MA48", "MA50", "MA52", "MA55", "MA56", "MA58", "MA60", "MA62", "MA64", "MA66", "MA68", "MA70", "11", "25", "33", "41", "43", "45", "47", "49", "51", "53", "55", "57", "59", "63", "65", "67", "69", "73", "75", "77", "79", "81", "83", "85", "87", "91", "93", "95", "97", "99", "101", "103", "105", "108", "111", "115", "117", "119", "173", "174", "175", "177", "178", "180", "181", "186", "187", "188", "189", "190", "191", "192", "194", "195", "197", "198", "199", "200", "203", "205", "207", "208", "Al0003", "Al0004", "Al0005", "Al0006", "Al0008", "Al0010", "Al0011", "Al0013", "Al0014", "Al0015", "Al0016", "Al0017", "Al0019", "Al0020", "Al0021", "Al0022", "Al0023", "Al0024", "Al0026", "Al0037", "Al0038", "Al0039", "Al0042", "Al0043", "Al0044", "Al0045", "Al0047", "Al0048", "Al0049", "Al0051", "Al0052", "Al0054", "Al0056", "Al0057", "Al0059", "Al0060", "Al0061", "Al0062", "Al0102", "Al0104", "Al0106", "Al0108", "Al0111", "Al0112", "Al0113", "Al0114", "Al0115", "Al0117", "Al0118", "Al0119", "Al0131", "Al0132", "Al0133", "Al0134", "Al0135", "Al0136", "Al0158", "Al0159", "Al0161", "Al0162", "Al0163", "Al0164", "Al0165", "Al0166", "Al0168", "Al0169", "Al0170", "Al0171", "Al0173", "Al0174", "Al0175", "Al0176", "Al0178", "Al0179", "Al0181", "Al0182", "Al0184", "Al0186", "Al0187", "Al0188", "Al0189", "Al0191", "Al0192", "Al0208", "Al0500", "Al0501", "Al0376", "Al0377", "Al0392", "Al0393", "Al0400", "Al0402", "Al0403", "Al0404", "Al0405", "Al0407", "Al0408", "Al0409", "Al0410", "Al0412", "Al0414", "Al0415", "Al0419", "Al0420", "Al0421", "Al0422", "Al0423", "Al0424", "Al0425", "Al0426", "Al0432", "Al0433", "Al0435", "Al0436", "Al0437", "Al0438", "Al0439", "Al0440", "Al0441", "Al0442", "Al0443", "Al0451", "Al0452", "Al0454", "Al0455", "Al0456", "Al0457", "Al0458", "Al0459", "Al0460", "Al0461", "Al0462", "Al0463", "Al0464", "Al0465", "Al0469", "Al0471", "Al0473", "Al0474", "Al0475", "Al0476", "Al0483", "Al0621", "Al0622", "Al0623", "Al0624", "Al0625", "Al0626", "Al0627", "Al0628", "Al0629", "Al0630", "Al0631", "Al0632", "Al0633", "Al0635", "Al0636", "Al0637", "Al0638", "Al0639", "Al0645", "Al0646", "Al0647", "Al0658", "Al0691", "Al0693", "Al0694", "Al0696", "Al0697", "Al0699", "Al0701", "Al0703", "Al0704", "Al0705", "Al0706", "Al0709", "Al0715", "Al0716", "Al0720", "Al0351", "Al0352", "Al0353", "Al0354", "Al0355", "Al0356", "Al0357", "Al0358", "Al0359", "Al0360", "Al0361", "Al0362", "Al0363", "Al0364", "Al0365", "Al0366", "Al0367", "Al0368", "Al0369", "Al0370", "Al0371", "Al0258", "Al0259", "Al0260", "Al0261", "Al0262", "Al0263", "Al0264", "Al0265", "Al0266", "Al0267", "Al0268", "Al0269", "Al0271", "Al0272", "Al0273", "Al0274", "Al0275", "Al0276", "Al0277", "Al0283", "Al0284", "Al0286", "Al0287", "Al0293", "Al0294", "Al0295", "Al0309", "Al0310", "Al0311", "Al0313", "Al0315", "Al0317", "Al0318", "Al0319", "Al0320", "Al0745", "Al0746", "Al0747", "Al0748", "Al0749", "Al0750", "Al0751", "Al0752", "Al0753", "Al0754", "Al0755", "Al0756", "Al0757", "Al0758", "Al0759", "Al0760", "Al0761", "Al0762", "Al0763", "Al0764", "Al0765", "Al0773", "Al0774", "Al0805", "IndA", "IndC", "IndD", "IndE", "IndF", "IndH", "IndI", "IndJ", "61", "113", "183", "193", "196", "202", "204", "206", "Al0002", "Al0009", "Al0101", "Al0103", "Al0105", "Al0137", "Al0138", "Al0139", "Al0140", "Al0141", "Al0142", "Al0144", "Al0145", "Al0146", "Al0148", "Al0149", "Al0150", "Al0152", "Al0153", "Al0154", "Al0155", "Al0156", "Al0177", "Al0180", "Al0183", "Al0185", "Al0210", "Al0214", "Al0220", "Al0224", "Al0230", "Al0232", "Al0236", "Al0502", "Al0416", "Al0417", "Al0418", "Al0444", "Al0446", "Al0447", "Al0448", "Al0279", "Al0282", "Al0285", "Al0288", "Al0289", "Al0290", "Al0292", "121", "122", "123", "125", "126", "127", "128", "129", "130", "131", "132", "133", "134", "135", "136", "137", "138", "139", "140", "141", "142", "143", "144", "145", "146", "147", "148", "149", "150", "151", "152", "153", "155", "156", "157", "158", "159", "161", "162", "163", "164", "165", "166", "167", "168", "170", "171", "172", "179", "184", "185", "274", "275", "277", "279", "281", "283", "285", "287", "289", "291", "353", "415", "417", "419", "421", "423", "425", "427", "429", "431", "433", "435", "437", "439", "441", "503", "505", "507", "509", "568", "570", "572", "574", "576", "578", "580", "582", "584", "586", "588", "590", "592", "594", "596", "598", "600", "602", "604", "718", "720", "722", "726", "728", "730", "732", "734", "736", "830", "832", "834", "836", "838", "840", "842", "844", "845", "847", "849", "851", "853", "855", "857", "1000", "1002", "1004", "1006", "1008", "1010", "1012", "1014", "1016", "1018", "1020", "1022", "1024", "1026", "1028", "1030", "1032", "1034", "1036", "1041", "1043", "1045", "1047", "1049", "1051", "1053", "1055", "1059", "1061", "1063", "1065", "1067", "1069", "1071", "1073", "1075", "1077", "1079", "1081", "1085", "1087", "1089", "1093", "1095")

pdf("Fractional_abundance_filtered_reordered.pdf", width=30, height=14)
p <- plot_bar(symbiontstrans, x= "Sample", fill = "Symbiont") + facet_grid(. ~ Host, scales = "free", space = "free") + geom_bar(aes(color=Symbiont, fill=Symbiont), stat='identity', position='stack') + scale_color_manual(values = col) + scale_fill_manual(values = col) + ylab("Fractional abundance") + theme_bw() + theme(axis.title = element_text(size=15, face="bold")) + theme(axis.text.x = element_blank(), axis.text.y = element_text(size=13)) + theme(legend.text = element_text(size = 13)) + theme(legend.title = element_text(size = 15, face="bold")) + theme(strip.text.x = element_blank()) + theme(legend.key = element_rect(size = 0.4)) + theme(legend.position="bottom")
p$data$Sample <- factor(p$data$Sample, levels = desired_order)
print(p)
dev.off()

colors <- c("pink1", "dodgerblue1", "goldenrod1", "red1", "maroon4", "cyan4", "midnightblue")
colors2 <- c("orangered3", "royalblue", "orange", "plum4", "darkseagreen", "darkgreen", "steelblue1", "gold", "lightsteelblue1", "plum3")

# PCoA ordination plots
PCoA1 <- ordinate(symbiontstrans, "PCoA", "unifrac", weighted = TRUE)
p1 <- plot_ordination(symbiontstrans, PCoA1)
pdf("PCoA_Unifrac_TSS.pdf")
p1 + scale_fill_manual(values = colors, name="Host") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Host)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15), legend.text = element_text(face = "italic"))
dev.off()

PCoA2 <- ordinate(hellinger.symbionts, "PCoA", "unifrac", weighted = TRUE)
p2 <- plot_ordination(hellinger.symbionts, PCoA2)
pdf("PCoA_Unifrac_HEL.pdf")
p2 + scale_fill_manual(values = colors, name="Host") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Host)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15), legend.text = element_text(face = "italic"))
dev.off()

# Alpha diversity
pdf("Alpha_diversity.pdf")
plot_richness(symbiontstrans, x="Host", measures=c("Shannon", "Simpson")) + scale_fill_manual(values = colors2, name="Geographic region") + theme_bw() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Basin)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(axis.text.x = element_text(face = "italic"))
dev.off()

# PCoA/NMDS plots within species
colorsb <- c("royalblue", "orange", "plum4", "steelblue1", "lightsteelblue1")
colorsb2 <- c("royalblue4", "cornflowerblue")
colorsb3 <- c("palevioletred4", "palevioletred1", "slateblue4", "lightpink", "slateblue1")

distb <- UniFrac(bouchetitrans, weighted = TRUE)
ordb <- metaMDS(distb, k=2, trymax=1000)
scoresb <- as.data.frame(scores(ordb, display = "sites"))
scoresb <- cbind(scoresb, Basin = sample_data(bouchetitrans)$Basin)
pdf("NMDS_Unifrac_boucheti_TSS.pdf")
pb <- ggplot(scoresb, aes(x=NMDS1, y=NMDS2, color=Basin))
pb + scale_fill_manual(values = colorsb, name="Geographic region") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Basin)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

ELSCb = subset_samples(bouchetitrans, Basin == "ELSC")
distb2 <- UniFrac(ELSCb, weighted = TRUE)
ordb2 <- metaMDS(distb2, k=2, trymax=1000)
scoresb2 <- as.data.frame(scores(ordb2, display = "sites"))
scoresb2 <- cbind(scoresb2, Vent = sample_data(ELSCb)$Vent)
pdf("NMDS_Unifrac_boucheti_ELSC_TSS.pdf")
pb2 <- ggplot(scoresb2, aes(x=NMDS1, y=NMDS2, color=Vent))
pb2 + scale_fill_manual(values = colorsb2, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

MBb = subset_samples(bouchetitrans, Basin == "Manus Basin")
distb3 <- UniFrac(MBb, weighted = TRUE)
ordb3 <- metaMDS(distb3, k=2, trymax=1000)
scoresb3 <- as.data.frame(scores(ordb3, display = "sites"))
scoresb3 <- cbind(scoresb3, Vent = sample_data(MBb)$Vent)
pdf("NMDS_Unifrac_boucheti_MB_TSS.pdf")
pb3 <- ggplot(scoresb3, aes(x=NMDS1, y=NMDS2, color=Vent))
pb3 + scale_fill_manual(values = colorsb3, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

distb <- UniFrac(hellinger.boucheti, weighted = TRUE)
ordb <- metaMDS(distb, k=2, trymax=1000)
scoresb <- as.data.frame(scores(ordb, display = "sites"))
scoresb <- cbind(scoresb, Basin = sample_data(hellinger.boucheti)$Basin)
pdf("NMDS_Unifrac_boucheti_HEL.pdf")
pb <- ggplot(scoresb, aes(x=NMDS1, y=NMDS2, color=Basin))
pb + scale_fill_manual(values = colorsb, name="Geographic region") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Basin)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

ELSCb = subset_samples(hellinger.boucheti, Basin == "ELSC")
distb2 <- UniFrac(ELSCb, weighted = TRUE)
ordb2 <- metaMDS(distb2, k=2, trymax=1000)
scoresb2 <- as.data.frame(scores(ordb2, display = "sites"))
scoresb2 <- cbind(scoresb2, Vent = sample_data(ELSCb)$Vent)
pdf("NMDS_Unifrac_boucheti_ELSC_HEL.pdf")
pb2 <- ggplot(scoresb2, aes(x=NMDS1, y=NMDS2, color=Vent))
pb2 + scale_fill_manual(values = colorsb2, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

MBb = subset_samples(hellinger.boucheti, Basin == "Manus Basin")
distb3 <- UniFrac(MBb, weighted = TRUE)
ordb3 <- metaMDS(distb3, k=2, trymax=1000)
scoresb3 <- as.data.frame(scores(ordb3, display = "sites"))
scoresb3 <- cbind(scoresb3, Vent = sample_data(MBb)$Vent)
pdf("NMDS_Unifrac_boucheti_MB_HEL.pdf")
pb3 <- ggplot(scoresb3, aes(x=NMDS1, y=NMDS2, color=Vent))
pb3 + scale_fill_manual(values = colorsb3, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()



colorsk <- c("royalblue", "orange", "plum4", "steelblue1", "gold", "plum3")
colorsk2 <- c("royalblue4", "cornflowerblue", "slategray1")
colorsk3 <- c("slateblue4", "plum", "slateblue1")

distk <- UniFrac(kojimaitrans, weighted = TRUE)
ordk <- metaMDS(distk, k=2, trymax=1000)
scoresk <- as.data.frame(scores(ordk, display = "sites"))
scoresk <- cbind(scoresk, Basin = sample_data(kojimaitrans)$Basin)
pdf("NMDS_Unifrac_kojimai_TSS.pdf")
pk <- ggplot(scoresk, aes(x=NMDS1, y=NMDS2, color=Basin))
pk + scale_fill_manual(values = colorsk, name="Geographic region") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Basin)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

ELSCk = subset_samples(kojimaitrans, Basin == "ELSC")
distk2 <- UniFrac(ELSCk, weighted = TRUE)
ordk2 <- metaMDS(distk2, k=2, trymax=1000)
scoresk2 <- as.data.frame(scores(ordk2, display = "sites"))
scoresk2 <- cbind(scoresk2, Vent = sample_data(ELSCk)$Vent)
pdf("NMDS_Unifrac_kojimai_ELSC_TSS.pdf")
pk2 <- ggplot(scoresk2, aes(x=NMDS1, y=NMDS2, color=Vent))
pk2 + scale_fill_manual(values = colorsk2, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pcoak <- cmdscale(distk2, add=TRUE, eig=TRUE)
Axis.1 <- pcoak$points[,1]
Axis.2 <- pcoak$points[,2]
eigs <- eigenvals(pcoak)
var1 <- round(eigs[1]/sum(eigs), digits=4)*100
var2 <- round(eigs[2]/sum(eigs), digits=4)*100
x_axis <- paste("Axis.1"," ","[",var1,"%]",sep="")
y_axis <- paste("Axis.2"," ","[",var2,"%]",sep="")
data <- cbind(as.data.frame(Axis.1), as.data.frame(Axis.2), Vent = sample_data(ELSCk)$Vent)
pdf("PCoA_Unifrac_kojimai_ELSC_TSS.pdf")
pk2b <- ggplot(data, aes(Axis.1, Axis.2, color=Vent))
pk2b + scale_fill_manual(values = colorsk2, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15)) + labs(x = x_axis, y = y_axis)
dev.off()

MBk = subset_samples(kojimaitrans, Basin == "Manus Basin")
distk3 <- UniFrac(MBk, weighted = TRUE)
ordk3 <- metaMDS(distk3, k=2, trymax=1000)
scoresk3 <- as.data.frame(scores(ordk3, display = "sites"))
scoresk3 <- cbind(scoresk3, Vent = sample_data(MBk)$Vent)
pdf("NMDS_Unifrac_kojimai_MB_TSS.pdf")
pk3 <- ggplot(scoresk3, aes(x=NMDS1, y=NMDS2, color=Vent))
pk3 + scale_fill_manual(values = colorsk3, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pcoak2 <- cmdscale(distk3, add=TRUE, eig=TRUE)
Axis.1 <- pcoak2$points[,1]
Axis.2 <- pcoak2$points[,2]
eigs <- eigenvals(pcoak2)
var1 <- round(eigs[1]/sum(eigs), digits=4)*100
var2 <- round(eigs[2]/sum(eigs), digits=4)*100
x_axis <- paste("Axis.1"," ","[",var1,"%]",sep="")
y_axis <- paste("Axis.2"," ","[",var2,"%]",sep="")
data <- cbind(as.data.frame(Axis.1), as.data.frame(Axis.2), Vent = sample_data(MBk)$Vent)
pdf("PCoA_Unifrac_kojimai_MB_TSS.pdf")
pk3b <- ggplot(data, aes(Axis.1, Axis.2, color=Vent))
pk3b + scale_fill_manual(values = colorsk3, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15)) + labs(x = x_axis, y = y_axis)
dev.off()

distk <- UniFrac(hellinger.kojimai, weighted = TRUE)
ordk <- metaMDS(distk, k=2, trymax=1000)
scoresk <- as.data.frame(scores(ordk, display = "sites"))
scoresk <- cbind(scoresk, Basin = sample_data(hellinger.kojimai)$Basin)
pdf("NMDS_Unifrac_kojimai_HEL.pdf")
pk <- ggplot(scoresk, aes(x=NMDS1, y=NMDS2, color=Basin))
pk + scale_fill_manual(values = colorsk, name="Geographic region") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Basin)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

ELSCk = subset_samples(hellinger.kojimai, Basin == "ELSC")
distk2 <- UniFrac(ELSCk, weighted = TRUE)
ordk2 <- metaMDS(distk2, k=2, trymax=1000)
scoresk2 <- as.data.frame(scores(ordk2, display = "sites"))
scoresk2 <- cbind(scoresk2, Vent = sample_data(ELSCk)$Vent)
pdf("NMDS_Unifrac_kojimai_ELSC_HEL.pdf")
pk2 <- ggplot(scoresk2, aes(x=NMDS1, y=NMDS2, color=Vent))
pk2 + scale_fill_manual(values = colorsk2, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pcoak <- cmdscale(distk2, add=TRUE, eig=TRUE)
Axis.1 <- pcoak$points[,1]
Axis.2 <- pcoak$points[,2]
eigs <- eigenvals(pcoak)
var1 <- round(eigs[1]/sum(eigs), digits=4)*100
var2 <- round(eigs[2]/sum(eigs), digits=4)*100
x_axis <- paste("Axis.1"," ","[",var1,"%]",sep="")
y_axis <- paste("Axis.2"," ","[",var2,"%]",sep="")
data <- cbind(as.data.frame(Axis.1), as.data.frame(Axis.2), Vent = sample_data(ELSCk)$Vent)
pdf("PCoA_Unifrac_kojimai_ELSC_HEL.pdf")
pk2b <- ggplot(data, aes(Axis.1, Axis.2, color=Vent))
pk2b + scale_fill_manual(values = colorsk2, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15)) + labs(x = x_axis, y = y_axis)
dev.off()

MBk = subset_samples(hellinger.kojimai, Basin == "Manus Basin")
distk3 <- UniFrac(MBk, weighted = TRUE)
ordk3 <- metaMDS(distk3, k=2, trymax=1000)
scoresk3 <- as.data.frame(scores(ordk3, display = "sites"))
scoresk3 <- cbind(scoresk3, Vent = sample_data(MBk)$Vent)
pdf("NMDS_Unifrac_kojimai_MB_HEL.pdf")
pk3 <- ggplot(scoresk3, aes(x=NMDS1, y=NMDS2, color=Vent))
pk3 + scale_fill_manual(values = colorsk3, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pcoak2 <- cmdscale(distk3, add=TRUE, eig=TRUE)
Axis.1 <- pcoak2$points[,1]
Axis.2 <- pcoak2$points[,2]
eigs <- eigenvals(pcoak2)
var1 <- round(eigs[1]/sum(eigs), digits=4)*100
var2 <- round(eigs[2]/sum(eigs), digits=4)*100
x_axis <- paste("Axis.1"," ","[",var1,"%]",sep="")
y_axis <- paste("Axis.2"," ","[",var2,"%]",sep="")
data <- cbind(as.data.frame(Axis.1), as.data.frame(Axis.2), Vent = sample_data(MBk)$Vent)
pdf("PCoA_Unifrac_kojimai_MB_HEL.pdf")
pk3b <- ggplot(data, aes(Axis.1, Axis.2, color=Vent))
pk3b + scale_fill_manual(values = colorsk3, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15)) + labs(x = x_axis, y = y_axis)
dev.off()


colorss <- c("royalblue", "orange", "gold")
colorss2 <- c("royalblue4", "cornflowerblue", "slategray1")

dists <- UniFrac(strummeritrans, weighted = TRUE)
ords <- metaMDS(dists, k=2, trymax=1000)
scoress <- as.data.frame(scores(ords, display = "sites"))
scoress <- cbind(scoress, Basin = sample_data(strummeritrans)$Basin)
pdf("NMDS_Unifrac_strummeri_TSS.pdf")
ps <- ggplot(scoress, aes(x=NMDS1, y=NMDS2, color=Basin))
ps + scale_fill_manual(values = colorss, name="Geographic region") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Basin)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

PCoAs <- ordinate(strummeritrans, "PCoA", "unifrac", weighted = TRUE)
psb <- plot_ordination(strummeritrans, PCoAs)
pdf("PCoA_Unifrac_strummeri_TSS.pdf")
psb + scale_fill_manual(values = colorss, name="Geographic region") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Basin)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

ELSCs = subset_samples(strummeritrans, Basin == "ELSC")
dists2 <- UniFrac(ELSCs, weighted = TRUE)
ords2 <- metaMDS(dists2, k=2, trymax=1000)
scoress2 <- as.data.frame(scores(ords2, display = "sites"))
scoress2 <- cbind(scoress2, Vent = sample_data(ELSCs)$Vent)
pdf("NMDS_Unifrac_strummeri_ELSC_TSS.pdf")
ps2 <- ggplot(scoress2, aes(x=NMDS1, y=NMDS2, color=Vent))
ps2 + scale_fill_manual(values = colorss2, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

dists <- UniFrac(hellinger.strummeri, weighted = TRUE)
ords <- metaMDS(dists, k=2, trymax=1000)
scoress <- as.data.frame(scores(ords, display = "sites"))
scoress <- cbind(scoress, Basin = sample_data(hellinger.strummeri)$Basin)
pdf("NMDS_Unifrac_strummeri_HEL.pdf")
ps <- ggplot(scoress, aes(x=NMDS1, y=NMDS2, color=Basin))
ps + scale_fill_manual(values = colorss, name="Geographic region") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Basin)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

PCoAs2 <- ordinate(hellinger.strummeri, "PCoA", "unifrac", weighted = TRUE)
psc <- plot_ordination(hellinger.strummeri, PCoAs2)
pdf("PCoA_Unifrac_strummeri_HEL.pdf")
psc + scale_fill_manual(values = colorss, name="Geographic region") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Basin)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

ELSCs = subset_samples(hellinger.strummeri, Basin == "ELSC")
dists2 <- UniFrac(ELSCs, weighted = TRUE)
ords2 <- metaMDS(dists2, k=2, trymax=1000)
scoress2 <- as.data.frame(scores(ords2, display = "sites"))
scoress2 <- cbind(scoress2, Vent = sample_data(ELSCs)$Vent)
pdf("NMDS_Unifrac_strummeri_ELSC_HEL.pdf")
ps2 <- ggplot(scoress2, aes(x=NMDS1, y=NMDS2, color=Vent))
ps2 + scale_fill_manual(values = colorss2, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()



colorsh <- c("darkslategray", "darkseagreen3", "khaki", "aquamarine4", "lightgoldenrodyellow")

disth <- UniFrac(hessleritrans, weighted = TRUE)
ordh <- metaMDS(disth, k=2, trymax=1000)
scoresh <- as.data.frame(scores(ordh, display = "sites"))
scoresh <- cbind(scoresh, Vent = sample_data(hessleritrans)$Vent)
pdf("NMDS_Unifrac_hessleri_TSS.pdf")
ph <- ggplot(scoresh, aes(x=NMDS1, y=NMDS2, color=Vent))
ph + scale_fill_manual(values = colorsh, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

disth <- UniFrac(hellinger.hessleri, weighted = TRUE)
ordh <- metaMDS(disth, k=2, trymax=1000)
scoresh <- as.data.frame(scores(ordh, display = "sites"))
scoresh <- cbind(scoresh, Vent = sample_data(hellinger.hessleri)$Vent)
pdf("NMDS_Unifrac_hessleri_HEL.pdf")
ph <- ggplot(scoresh, aes(x=NMDS1, y=NMDS2, color=Vent))
ph + scale_fill_manual(values = colorsh, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()


colorsn <- c("royalblue", "lightsteelblue1")
colorsn2 <- c("royalblue4", "steelblue", "slategray1")

distn <- UniFrac(nautileitrans, weighted = TRUE)
ordn <- metaMDS(distn, k=2, trymax=1000)
scoresn <- as.data.frame(scores(ordn, display = "sites"))
scoresn <- cbind(scoresn, Basin = sample_data(nautileitrans)$Basin)
pdf("NMDS_Unifrac_nautilei_TSS.pdf")
pn <- ggplot(scoresn, aes(x=NMDS1, y=NMDS2, color=Basin))
pn + scale_fill_manual(values = colorsn, name="Geographic region") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Basin)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

ELSCn = subset_samples(nautileitrans, Basin == "ELSC")
distn2 <- UniFrac(ELSCn, weighted = TRUE)
ordn2 <- metaMDS(distn2, k=2, trymax=1000)
scoresn2 <- as.data.frame(scores(ordn2, display = "sites"))
scoresn2 <- cbind(scoresn2, Vent = sample_data(ELSCn)$Vent)
pdf("NMDS_UniFrac_nautilei_ELSC_TSS.pdf")
pn2 <- ggplot(scoresn2, aes(x=NMDS1, y=NMDS2, color=Vent))
pn2 + scale_fill_manual(values = colorsn2, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

distn <- UniFrac(hellinger.nautilei, weighted = TRUE)
ordn <- metaMDS(distn, k=2, trymax=1000)
scoresn <- as.data.frame(scores(ordn, display = "sites"))
scoresn <- cbind(scoresn, Basin = sample_data(hellinger.nautilei)$Basin)
pdf("NMDS_Unifrac_nautilei_HEL.pdf")
pn <- ggplot(scoresn, aes(x=NMDS1, y=NMDS2, color=Basin))
pn + scale_fill_manual(values = colorsn, name="Geographic region") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Basin)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

ELSCn = subset_samples(hellinger.nautilei, Basin == "ELSC")
distn2 <- UniFrac(ELSCn, weighted = TRUE)
ordn2 <- metaMDS(distn2, k=2, trymax=1000)
scoresn2 <- as.data.frame(scores(ordn2, display = "sites"))
scoresn2 <- cbind(scoresn2, Vent = sample_data(ELSCn)$Vent)
pdf("NMDS_UniFrac_nautilei_ELSC_HEL.pdf")
pn2 <- ggplot(scoresn2, aes(x=NMDS1, y=NMDS2, color=Vent))
pn2 + scale_fill_manual(values = colorsn2, name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()


# PERMANOVA and LDM
sampledf <- data.frame(sample_data(symbionts))
otudf <- t(otu_table(symbionts))
tree = phy_tree(symbionts)

# Model 1
fit <- ldm(otudf|(Preservation+Extraction+Sequencing) ~ Host + Basin, data=sampledf, n.perm.max=10000, seed=12345, scale.otu.table=TRUE, dist.method="wt-unifrac", tree=tree, n.cores=4)

fit$F.global.freq
fit$F.global.tran
fit$p.global.freq
fit$p.global.tran
fit$p.global.omni
fit$VE.df.submodels

pdf("LDM_screeplot.pdf")
par(mfrow = c(1,2))
scree.freq <- c(fit$VE.global.freq.submodels/fit$VE.df.submodels,fit$VE.global.freq.residuals)
color <- c("red2", "midnightblue", rep("gray50", length(scree.freq)-2))
plot(scree.freq/sum(scree.freq), main="Frequency Scale", xlab="Component", ylab="Proportion of total sum of squares", col=color, pch=16)
scree.tran <- c(fit$VE.global.tran.submodels/fit$VE.df.submodels,fit$VE.global.tran.residuals)
color <- c("red2", "midnightblue", rep("gray50", length(scree.tran)-2))
plot(scree.tran/sum(scree.tran), main="Arcsin-Root Scale", xlab="Component", ylab="", col=color, pch=16)
dev.off()

round(scree.freq/sum(scree.freq)*100, 2)[1]
round(scree.freq/sum(scree.freq)*100, 2)[2]

perm <- permanovaFL(otudf|(Preservation+Extraction+Sequencing) ~ Host + Basin, data=sampledf, n.perm.max=1000, seed=12345, scale.otu.table=TRUE, dist.method="wt-unifrac", tree=tree, n.cores=4)

perm$F.statistics
perm$p.permanova

# Subset data
subset <- subset_samples(symbionts, Basin == "ELSC")
subset = subset_samples(subset, Host == "A. boucheti" | Host == "A. kojimai" | Host == "A. strummeri")
subset = filter_taxa(subset, function(x) sum(x) > 0, TRUE)
subsettrans <- transform_sample_counts(subset, function(x) x/sum(x))

sampledf2 <- data.frame(sample_data(subset))
otudf2 <- t(otu_table(subset))
tree2 <- phy_tree(subset)

# Model 2
fit2 <- ldm(otudf2|(Preservation+Extraction) ~ Vent + Host, data=sampledf2, n.perm.max=10000, seed=12345, scale.otu.table=TRUE, dist.method="wt-unifrac", tree=tree2, n.cores=4)

fit2$F.global.freq
fit2$F.global.tran
fit2$p.global.freq
fit2$p.global.tran
fit2$p.global.omni
fit2$VE.df.submodels

pdf("LDM_screeplot2.pdf")
par(mfrow = c(1,2))
scree.freq2 <- c(fit2$VE.global.freq.submodels/fit2$VE.df.submodels,fit2$VE.global.freq.residuals)
color <- c("red2", "midnightblue", rep("gray50", length(scree.freq2)-2))
plot(scree.freq2/sum(scree.freq2), main="Frequency Scale", xlab="Component", ylab="Proportion of total sum of squares", col=color, pch=16)
scree.tran <- c(fit2$VE.global.tran.submodels/fit2$VE.df.submodels,fit2$VE.global.tran.residuals)
color <- c("red2", "midnightblue", rep("gray50", length(scree.tran)-2))
plot(scree.tran/sum(scree.tran), main="Arcsin-Root Scale", xlab="Component", ylab="", col=color, pch=16)
dev.off()

round(scree.freq2/sum(scree.freq2)*100, 2)[1]
round(scree.freq2/sum(scree.freq2)*100, 2)[2]

perm2 <- permanovaFL(otudf2|(Preservation+Extraction) ~ Vent + Host, data=sampledf2, n.perm.max=1000, seed=12345, scale.otu.table=TRUE, dist.method="wt-unifrac", tree=tree2, n.cores=4)

perm2$F.statistics
perm2$p.permanova

# Model 3
fit3 <- ldm(otudf2 ~ Preservation*Extraction + Vent + Host, data=sampledf2, n.perm.max=10000, seed=12345, scale.otu.table=TRUE, dist.method="wt-unifrac", tree=tree2, n.cores=4)

fit3$F.global.freq
fit3$F.global.tran
fit3$p.global.freq
fit3$p.global.tran
fit3$p.global.omni
fit3$VE.df.submodels

pdf("LDM_screeplot3.pdf")
par(mfrow = c(1,2))
scree.freq3 <- c(fit3$VE.global.freq.submodels/fit3$VE.df.submodels,fit3$VE.global.freq.residuals)
color <- c("red2", "midnightblue", "lightseagreen", rep("gray50", length(scree.freq3)-3))
plot(scree.freq3/sum(scree.freq3), main="Frequency Scale", xlab="Component", ylab="Proportion of total sum of squares", col=color, pch=16)
scree.tran <- c(fit3$VE.global.tran.submodels/fit3$VE.df.submodels,fit3$VE.global.tran.residuals)
color <- c("red2", "midnightblue", "lightseagreen", rep("gray50", length(scree.tran)-3))
plot(scree.tran/sum(scree.tran), main="Arcsin-Root Scale", xlab="Component", ylab="", col=color, pch=16)
dev.off()

round(scree.freq3/sum(scree.freq3)*100, 2)[1]
round(scree.freq3/sum(scree.freq3)*100, 2)[2]
round(scree.freq3/sum(scree.freq3)*100, 2)[3]

perm3 <- permanovaFL(otudf2 ~ Preservation*Extraction + Vent + Host, data=sampledf2, n.perm.max=1000, seed=12345, scale.otu.table=TRUE, dist.method="wt-unifrac", tree=tree2, n.cores=4)

perm3$F.statistics
perm3$p.permanova


distsub <- UniFrac(subsettrans, weighted = TRUE)
ordsub <- metaMDS(distsub, k=2, trymax=1000)
scoressub <- as.data.frame(scores(ordsub, display = "sites"))
scoressub <- cbind(scoressub, Host = sample_data(subset)$Host, Vent = sample_data(subset)$Vent)
pdf("NMDS_Unifrac_subset_TSS.pdf")
psub <- ggplot(scoressub, aes(x=NMDS1, y=NMDS2, color=Host, shape=Vent))
psub + scale_fill_manual(values = c("dodgerblue1", "red1", "cyan4"), name="Host") + scale_shape_manual(values = c(21, 22, 24), name="Vent") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Host), shape = factor(Vent)), color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

scoressubm <- cbind(scoressub, Preservation = sample_data(subset)$Preservation, Extraction = sample_data(subset)$Extraction)
pdf("NMDS_Unifrac_Methods_subset_TSS.pdf")
psub <- ggplot(scoressubm, aes(x=NMDS1, y=NMDS2, color=Extraction, shape=Preservation))
psub + scale_fill_manual(values = c("lightseagreen", "slateblue"), name="Extraction") + scale_shape_manual(values = c(21, 24), name="Preservation") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Extraction), shape = factor(Preservation)), color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

