#!/bin/bash
#SBATCH -J Oligotyping
#SBATCH -t 24:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=300g
#SBATCH -o Oligotyping.out
#SBATCH -e Oligotyping.err

export DISPLAY=:0.0
export MPLBACKEND="agg"

# After identifying the symbiont ASVs based on taxonomic assignment create a list of zOTU IDs and extract read information from the zMap
grep -w -f symbiont_zotus.txt zmap_symbionts.txt | sed -e "s/ .*//g" > symbiont_readIDs.txt
seqtk subseq symbionts_merged.fq symbiont_readIDs.txt > snails.fq
convert_fastaqual_fastq.py -f snails.fq -c fastq_to_fastaqual
sed -i "s/>/>Sample-/g" snails.fna
sed -i "s/[.]/_/g" snails.fna
sed -i "s/>/>Sample-/g" snails.qual
sed -i "s/[.]/_/g" snails.qual
o-get-sample-info-from-fasta snails.fna 
o-pad-with-gaps snails.fna -o snails_pads.fa
entropy-analysis snails_pads.fa --qual-scores-file snails.qual --no-display
# Entropy components need to be identified based on visual analysis of the resulting entropy graphs; these are the final ones selected for oligotyping
oligotype snails_pads.fa snails_pads.fa-ENTROPY -o snails -C 3,12,15,20,42,44,48,55,56,57,58,60,61,76,77,90,92,111,112,113,115,116,128,130,138,174,175,183,204,205,209,228,260,264,270,274,277,283,289,30
0,305,306,307,308,309,310,313,314,315,316,319,326,327,328,337,343,344,363,370,373 --qual-scores-dict snails.qual.cPickle --qual-stats-dict snails.qual.STATS.cPickle -q 20 -M 100 -A 50 -a 2.37 -s 5 --no-
display -N 24 --generate-sets

export LC_ALL=en_US.utf8
export LANG=en_US.utf8

source activate qiime2-2021.4

qiime tools import --input-path snails/OLIGO-REPRESENTATIVES.fasta --output-path oligo-seqs.qza --type 'FeatureData[Sequence]'

qiime feature-classifier classify-sklearn --i-classifier silva-132-99-515F-926R-classifier.qza --i-reads oligo-seqs.qza --o-classification taxonomy_oligos.qza --p-n-jobs 1

qiime alignment mafft --i-sequences oligo-seqs.qza --o-alignment aligned_oligos.qza
qiime alignment mask --i-alignment aligned_oligos.qza --o-masked-alignment masked-aligned_oligos.qza
qiime phylogeny iqtree-ultrafast-bootstrap --i-alignment masked-aligned_oligos.qza --o-tree unrooted-tree_oligos.qza --p-perturb-nni-strength 0.2 --p-stop-iter 200 --p-n-runs 5 --p-bootstrap-replicates 5000 --p-n-cores 1 --verbose
qiime phylogeny midpoint-root --i-tree unrooted-tree_oligos.qza --o-rooted-tree rooted-tree_oligos.qza

qiime tools export --input-path rooted-tree_oligos.qza --output-path oligos
qiime tools export --input-path taxonomy_oligos.qza --output-path oligos

