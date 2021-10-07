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
symbionts = filter_taxa(symbionts, function(x) sum(x) > 0, TRUE)
symbiontstrans <- transform_sample_counts(symbionts, function(x) x/sum(x))

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
plot(spec, col = "black", lwd=2, ci = 0, xlab = "# individuals", ylab = "# ASVs")
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

# PCoA ordination plots
PCoA <- ordinate(symbiontstrans, "PCoA", "unifrac", weighted = TRUE)
pbc1 <- plot_ordination(symbiontstrans, PCoA)
pdf("PCoA_Unifrac_filtered.pdf")
pbc1 + scale_fill_manual(values = colors, name="Host") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Host)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15), legend.text = element_text(face = "italic"))
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
plot_richness(symbiontstrans, x="Host", measures=c("Shannon", "Simpson")) + scale_fill_manual(values = colors2, name="Basin") + theme_bw() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Basin)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "italic"))
dev.off()
