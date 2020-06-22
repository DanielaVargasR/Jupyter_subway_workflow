# Amplicon Sequence Variants with DADA2 in QIIME2 using GreenGenes and also SILVA taxonomic database

# Used demultipled and pairend files
# ASV processed in qiime2-2019.7
# Built a manifest.txt file to link the files with the real names of the samples because (saved as .tsv)
# Eg.

sample-id	forward-absolute-filepath	 reverse-absolute-filepath
A42	LA2ME2SS100_S53_L001_R1_001.fastq.gz	LA2ME2SS100_S53_L001_R2_001.fastq.gz
A43	LA2ME2SS101_S54_L001_R1_001.fastq.gz	LA2ME2SS101_S54_L001_R2_001.fastq.gz
A44	LA2ME2SS102_S55_L001_R1_001.fastq.gz	LA2ME2SS102_S55_L001_R2_001.fastq.gz

# Note: input format: Phred33

conda activate qiime2-2019.7

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path pe-64-manifest.txt \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
  
# Generate a file to  vizualize
qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization demux.qzv  

# Vizualize
  qiime tools view demux.qzv

# Denoising and filtering with DADA2
# The truncation and trimming sizes according to the quality of the F and R reads.
# I removed the first 9 bases of read F because they show "bias" in the nucleotide composition eg 100% A or C or G or T

  qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --p-trunc-len-f 260 \
  --p-trunc-len-r 220 \
  --p-trim-left-f 9 \
  --p-trim-left-r 0 \
  --o-table denoise_table_dada2.qza \
  --o-representative-sequences denoise_rep_set_dada2.qza \
  --o-denoising-stats denoise_stat_dada2.qza
  
# Resume feature table   
  qiime feature-table summarize \
  --i-table denoise_table_dada2.qza \
  --o-visualization denoise_table_dada2.qzv \
  --m-sample-metadata-file metadata_R112_dada2.txt

# Rep set dada2
qiime feature-table tabulate-seqs \
  --i-data denoise_rep_set_dada2.qza \
  --o-visualization rep-seqs.qzv
  
# Vizualization
  qiime metadata tabulate \
  --m-input-file denoise_stat_dada2.qza \
  --o-visualization denoise_stat_dada2.qzv

##################################
#######GREENGENES#################
##################################

# Import GreenGenes pre trained
# Downloaded all gg 13.8 folder
  qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path 99_otus.fasta \
  --output-path 99_otus.qza

# Upload 99_otu_taxonomy sequences not aligned
  qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path 99_otu_taxonomy.txt \
  --output-path ref-taxonomy.qza
  
#Now we cut with primers V3-V4
#341F (CCTACGGGNGGCWGCAG) and 805R (GTGGACTACHVGGGTWTCTAAT)
qiime feature-classifier extract-reads \
  --i-sequences 99_otus.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GTGGACTACHVGGGTWTCTAAT \
  --o-reads ref-seqs.qza
  
# Train the classifier
# We can now train a Naive Bayes classifier as follows, using the reference reads and taxonomy that we just created.
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier classifier.qza

# Classification  
  qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads denoise_rep_set_dada2.qza \
  --o-classification taxonomy.qza 
  
# Rarificacion at 4907 seqs/samples keep 100% of samples
# Export
  mkdir exports
  
  qiime tools export \
--input-path denoise_table_dada2.qza \
--output-path exports/feature-table

biom convert \
-i exports/feature-table/feature-table.biom \
-o exports/feature-table.tsv \
--to-tsv

##################################
##############SILVA###############
##################################

#V3-V4 classifier: https://github.com/Jiung-Wen/q2-silva-V3V4classifier
  
  qiime feature-classifier classify-sklearn   \
  --i-classifier silva_132_99_v3v4.qza  \ 
  --i-reads denoise_rep_set_dada2.qza   \
  --p-pre-dispatch 6   \
  --p-reads-per-batch 100 \  
  --p-n-jobs 2   \
  --o-classification ref-taxonomy_silva_132_99_v3v4.qza


