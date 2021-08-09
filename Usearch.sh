#!/bin/bash
#SBATCH -J Usearch
#SBATCH -t 48:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=100g
#SBATCH -o Usearch.out
#SBATCH -e Usearch.err

for file in `cat filelist.txt`
do
  java -jar /gpfs/data/rbeinart/Software/trinityrnaseq-v2.11.0/trinity-plugins/Trimmomatic/trimmomatic-0.36.jar PE -threads 24 -phred33 Raw_reads/${file}*R1_001.fastq.gz Raw_reads/${file}*R2_001.fastq.gz ${file}_paired_R1.fq ${file}_singles_R1.fq ${file}_paired_R2.fq ${file}_singles_R
2.fq ILLUMINACLIP:16S_adaptors.fa:2:30:10
done

mkdir Meta-analysis
cd Meta-analysis

usearch11 -fastq_mergepairs ../*_paired_R1.fq -fastqout symbionts_merged.fq -fastq_minmergelen 300 -fastq_maxmergelen 400 -relabel @
usearch11 -fastq_filter symbionts_merged.fq -fastq_maxee_rate 0.001 -fastaout symbionts_filtered.fa -fastq_minlen 300 -fastq_truncqual 20
usearch11 -fastx_uniques symbionts_filtered.fa -fastaout symbionts_uniques.fa -sizeout -relabel Uniq
usearch11 -unoise3 symbionts_uniques.fa -zotus zotus_symbionts.fa -tabbedout unoise3.txt
usearch11 -otutab symbionts_merged.fq -zotus zotus_symbionts.fa -sample_delim . -otutabout zotu_symbionts.txt -mapout zmap_symbionts.txt -maxrejects 1000 -strand both

export LC_ALL=en_US.utf8
export LANG=en_US.utf8

source activate qiime2-2019.10

biom convert -i zotu_symbionts.txt -o zotu-table.biom --to-hdf5 --table-type="OTU table"

qiime tools import --input-path zotu-table.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path zotu-table.qza
qiime tools import --input-path zotus_symbionts.fa --output-path zotu_seqs.qza --type 'FeatureData[Sequence]'

qiime tools import --type 'FeatureData[Sequence]' --input-path silva_132_99_16S.fna --output-path silva_132_99_otus.qza
qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path silva_132_99_16S_consensus_taxonomy.txt --output-path silva_132_99-taxonomy.qza
qiime feature-classifier extract-reads --i-sequences silva_132_99_otus.qza --p-f-primer GTGYCAGCMGCCGCGGTAA --p-r-primer CCGYCAATTYMTTTRAGTTT --p-min-length 300 --p-max-length 500 --o-reads silva_132_99-ref-seqs.qza
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads silva_132_99-ref-seqs.qza --i-reference-taxonomy silva_132_99-taxonomy.qza --o-classifier silva-132-99-515F-926R-classifier.qza

qiime feature-classifier classify-sklearn --i-classifier silva-132-99-515F-926R-classifier.qza --i-reads zotu_seqs.qza --o-classification taxonomy.qza --p-n-jobs 1

qiime alignment mafft --i-sequences zotu_seqs.qza --o-alignment aligned_symbionts.qza
qiime alignment mask --i-alignment aligned_symbionts.qza --o-masked-alignment masked-aligned_symbionts.qza
qiime phylogeny fasttree --i-alignment masked-aligned_symbionts.qza --o-tree unrooted-tree.qza
qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

qiime tools export --input-path zotu-table.qza --output-path unfiltered
qiime tools export --input-path zotu_seqs.qza --output-path unfiltered
qiime tools export --input-path rooted-tree.qza --output-path unfiltered
qiime tools export --input-path taxonomy.qza --output-path unfiltered
