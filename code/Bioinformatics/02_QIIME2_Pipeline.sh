#Processing eDNA metabarcode data using QIIME2 and cutadapt, followed by DADA2 taxonomy assignments using custom fish reference databases
#The qiime tutorials are useful and found at https://docs.qiime2.org/2022.2/tutorials/overview/#useful-points-for-beginners 
#First activate QIIME if it hasn't been, can also reactivate qiime if you close the window 
conda activate qiime2-2023.5 &&
source tab-qiime #activate tab completion

#check currently active conda environment
conda info

# View plugin, used to view any .qzv file in html and export tsv and fasta files
qiime tools view /path_to_file/filename.qzv

#Now import our data using a 'manifest' file of all fastq file names
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path pe33-COImanifest \
--input-format PairedEndFastqManifestPhred33V2 \
--output-path TraskCOI-combined-demux.qza


#check out the data for visualization
qiime demux summarize \
  --i-data TraskCOI-combined-demux.qza \
  --o-visualization COI-demux-vis.qzv ##save tsv file of per-sample-fastq-counts.tsv for optional step below ##

  
 ## OPTIONAL: filter out samples with less than 100 reads (can set this to any number) ##
qiime demux filter-samples \
  --i-demux TraskCOI-combined-demux.qza \
  --m-metadata-file /path_to_output_folder/per-sample-fastq-counts.tsv \
  --p-where 'CAST([forward sequence count] AS INT) > 100' \
  --o-filtered-demux /path_to_output_folder/filename_greater100reads.qza

#Now trim primers
qiime cutadapt trim-paired \
--i-demultiplexed-sequences TraskCOI-combined-demux.qza \
--p-cores 40 \
--p-front-f GGWACWGGWTGAACWGTWTAYCCYCC \
--p-front-r TAIACYTCIGGRTGICCRAARAAYCA \
--p-error-rate 0.11 \
--p-discard-untrimmed \
--p-match-read-wildcards \
--p-match-adapter-wildcards \
--p-minimum-length 40 \
--o-trimmed-sequences COI-demux-trimmed.qza \
--output-dir  trimmed \
--verbose

#Hopefully if you put the primer sequences in properly, >90% of reads will have a primer trimmed which will show while cutadapt runs


#denoise using dada2 which infers ASVs 
#Note: using --p-n-threads = 0 will use all threads available 
### for 16S use --p-trunc-len-f 125 and --p-trunc-len-r 125; 12S use 116 and 108 ###
# can add --p-min-overlap 12 or some other number if need be
qiime dada2 denoise-paired \
--i-demultiplexed-seqs COI-demux-trimmed.qza \
--p-trunc-len-f  220 \
--p-trunc-len-r  220 \
--p-n-threads 0 \
--p-n-reads-learn 1000000 \
--p-pooling-method independent \
--output-dir trimmed/dada2out \
--verbose


#Generate summaries of denoising stats and feature table
qiime feature-table summarize \
  --i-table dada2out-test/table.qza \
  --o-visualization dada2out-test/table.qzv \
  --m-sample-metadata-file ../Trask-sample-metadata.tsv &&
qiime feature-table tabulate-seqs \
  --i-data dada2out-test/representative_sequences.qza \
  --o-visualization dada2out-test/rep-seqs.qzv &&
qiime metadata tabulate \
  --m-input-file dada2out-test/denoising_stats.qza \
  --o-visualization dada2out-test/denoising-stats.qzv
  
qiime tools view /path_to_output_folder/filename_rep_seqs.qzv  ## export the ASV fasta file from the view for input into FuzzyID2 and BLAST


 ### export results to biom formatted file
qiime tools export \
--input-path dada2out-test/table.qza \
--output-path dada2out-test/Trask_filtered_table_biom ##specifying a folder output here, this tool will automatically export a file called 'feature-table.biom' to this folder

### convert biom to tsv
biom convert -i dada2out-test/Trask_filtered_table_biom/feature-table.biom \
-o dada2out-test/ESI16S_filtered_table_biom/Trask_feature_table_export.tsv \
--to-tsv

### OPTIONAL filtering after exporting to tsv
## Remove rare ASV's by calculating if an ASV has a read number that is less than 0.1% of the total read number of that ASV across all samples. 
## This is summing across columns in the exported feature table, calculating 0.1% of that sum, and removing all instances where read numbers were less than that number.
 
 #Generate a phylogenetic tree from our data
 cd dada2out/
 qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences representative_sequences.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
  
  #now use the rooted tree to generate some biodiversity stats such as dissimilarity matrices and PCoA (Principal coordinates analysis)
  qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 1500 \
  --p-n-jobs-or-threads auto \
  --m-metadata-file ../../2021-sample-metadata.tsv \
  --output-dir 16S-core-metrics-results
 
###################################################
######TAXONOMY - Nick did not run this part but
 ##ran the RDP classifier separately in the Linux 
 ##terminal but below is code that could be used 
 ##within Qiime2
##############################################
## there are lots of ways to do taxonomy, including just blasting, or building a reference database using rescript (below), or using FuzzyID2 with a custom library
  #using rescript to train our classifier
  qiime rescript filter-taxa \
  --i-taxonomy fish-16S-ref-tax.qza \
  --m-ids-to-keep-file fish-16S-ref-seqs-keep.qza \
  --o-filtered-taxonomy fish-16S-ref-taxa-keep.qza
  
  qiime rescript evaluate-taxonomy \
 --i-taxonomies fish-16S-ref-taxa-keep.qza \
 --o-taxonomy-stats fish-16S-ref-tax-keep-eval.qzv
 
 qiime metadata tabulate \
 --m-input-file fish-16S-ref-taxa-keep.qza \
 --o-visualization fish-16S-ref-tax-keep.qzv &&
 qiime rescript evaluate-seqs \
 --i-sequences fish-16S-ref-seqs-keep.qza \
 --p-kmer-lengths 32 16 8 \
 --o-visualization fish-16S-ref-seqs-keep-eval.qzv
 
 #Build and evaluate classifier
 #here, the --o-classifier output is of type TaxonomicClassifier and the -o-observed-taxonomy is FeatureData[Taxonomy] (same as --i-taxonomy)
 qiime rescript evaluate-fit-classifier \
 --i-sequences fish-16S-ref-seqs-keep.qza \
 --i-taxonomy fish-16S-ref-taxa-keep.qza \
 --p-n-jobs -1 \
 --o-classifier ncbi-16S-fish-refseqs-classifier.qza \
 --o-evaluation ncbi-16S-fish-refseqs-classifier-evaluation.qzv \
 --o-observed-taxonomy ncbi-16S-fish-refseqs-predicted-taxonomy.qza \
 --output-dir 16S-Classifier \
 --verbose 
 
 qiime rescript evaluate-taxonomy \
 --i-taxonomies fish-16S-ref-taxa-keep.qza ncbi-16S-fish-refseqs-predicted-taxonomy.qza \
 --p-labels ref-taxonomy predicted-taxonomy \
 --o-taxonomy-stats 16S-ref-taxonomy-evaluation.qzv \
 --verbose
 
#Now back to qiime to do our taxonomy
  qiime feature-classifier classify-sklearn \
  --i-classifier ../../../../ReferenceData/ncbi-16S-fish-refseqs-classifier.qza \
  --i-reads representative_sequences.qza \
  --o-classification 16S-taxonomy.qza

qiime metadata tabulate \
  --m-input-file 16S-taxonomy.qza \
  --o-visualization 16S-taxonomy.qzv
  



