taskID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
arg = commandArgs(trailingOnly=TRUE)
cell_line = as.character(arg)

library(pacman)
p_load(ChIPseeker, TxDb.Hsapiens.UCSC.hg19.knownGene, BSgenome.Hsapiens.UCSC.hg19, clusterProfiler, rtracklayer,
GenomicRanges, ENCODExplorer, ggplot2, dplyr, reshape, heatmaps, abind)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
genome = BSgenome.Hsapiens.UCSC.hg19
ENCODE_promoters = readRDS("../data/ENCODE_promoters.rds") 
### extend promoter regions to 2.5kb
ENCODE_promoters_ext = lapply(ENCODE_promoters, function(x){ranges(x)+2250})
ENCODE_promoters_ext = ENCODE_promoters_ext[-5]
ENCODE_promoters_ext = lapply(seq_along(ENCODE_promoters_ext),function(x){ENCODE_promoters_ext[[x]]=GRanges(seqnames=seqnames(ENCODE_promoters[[x]]),
ranges=ENCODE_promoters_ext[[x]], seqlengths=seqlengths(genome))})
chr_set = paste("chr", c(1:22), sep = "")
for(i in 1:length(ENCODE_promoters_ext)){seqlevels(ENCODE_promoters_ext[[i]], pruning.mode="coarse")=chr_set}

histone_marks = c("H3K4me3", "H3K36me3", "H3K4me1", "H3K27me3", "H3K9ac", "H2AFZ", "H4K20me1", "H3K9me3", "H3K4me2", "H3K27ac", "H3K79me2", "H3K9me1")
cell_lines = c("K562", "HepG2", "A549", "GM12878")
names(ENCODE_promoters_ext) = cell_lines
saveRDS(ENCODE_promoters_ext, "../data/ENCODE_promoters_ext.rds")
##############################

### ENCODEexplorer package for data retrieval
data(encode_df, package = "ENCODExplorer")

.query_results = function(x){
   query_results = queryEncode(df=encode_df, organism = "Homo sapiens", biosample_name = cell_lines[i], target = histone_marks[x], assay = "ChIP-seq",
                               file_format = "bam", fixed = TRUE)
   query_results = query_results[assembly=="hg19"]
   return(query_results)
}

.query_results = function(x){
   query_results = queryEncode(df=encode_df, organism = "Homo sapiens", biosample_name = cell_lines[i], target = histone_marks[[x]], assay = "ChIP-seq", file_format = "bam", fixed = TRUE)
   query_results = query_results[assembly=="hg19"]
   return(query_results)
}

qr_list = list()
for(i in 1:length(cell_lines)){
   query_results = lapply(seq_along(histone_marks), .query_results)
   query_results = do.call(args=query_results, what=rbind)
   qr_list[[i]] = query_results
}
qr_list = do.call(args=qr_list, what=rbind)
# drop MCF-7 cell line
qr_list = qr_list[biosample_name != "MCF-7"]

## retrieve control files
.control_results = function(x){
   control_results = queryEncode(df=encode_df, organism = "Homo sapiens", biosample_name = cell_lines[x], target = "Control",
                                 assay = "ChIP-seq", file_format = "bam", fixed = TRUE)
   control_results = control_results[assembly=="hg19"]
   return(control_results)
}
control_list = lapply(seq_along(cell_lines), .control_results)
control_list = do.call(args=control_list, what=rbind)
# keep control files for our samples only
control_list = control_list[accession %in% qr_list$controls]

acc = c("ENCSR000AKY", "ENCSR000AME", "ENCSR000ASS", "ENCSR000AKJ")
control_list = control_list[accession %in% acc]
qr_list = rbind(qr_list, control_list)

################################################

p_load(GenomicAlignments)
# remove samples without index
accession_codes = qr_list[target != "H3K9me1"]$file_accession
accession_codes = accession_codes[!(accession_codes %in% c("ENCFF386DOZ", "ENCFF001EXQ", "ENCFF000AHU", "ENCFF000ASE",
                                 "ENCFF000ARA", "ENCFF000AHC", "ENCFF000BDP"))]
# construct a df of corresponding cell lines
.cl = function(x){
   cell_type = qr_list[file_accession==accession_codes[[x]]]$biosample_name
   return(cell_type)
}
.cl_1 = function(x){
   target = qr_list[file_accession==accession_codes[[x]]]$target
   return(target)
}
cell_type = as.character(lapply(seq_along(accession_codes), .cl))
target = as.character(lapply(seq_along(accession_codes), .cl_1))
cl = data.frame(cbind(cell_type, target), stringsAsFactors = FALSE)
cl["accession"] = accession_codes
cl$cell_type = as.character(cl$cell_type)
cl$target = as.character(cl$target)
saveRDS(cl, "../data/hmm_training/accession_df.rds")
hm = readRDS("../data/promoters_HM.rds")


### prepare data for chomatin state learning
ENCODE_promoters_ext = ENCODE_promoters_ext[[cell_line]]
cl_subset = subset(cl, cell_type==cell_line)

.binned_avg = function(my_promoters){
bins_promoters = tile(my_promoters, 50)
hm_subs = hm[[cl_subset$accession[taskID]]] # <----- branching step 
overlaps = lapply(bins_promoters, function(x){countOverlaps(x, hm_subs)})
overlaps = do.call(rbind, overlaps)
return(overlaps)
}
overlap_mat = .binned_avg(ENCODE_promoters_ext)
path = paste0("../data/hmm_training/matrices/", arg, "_", taskID, ".rds")
saveRDS(overlap_mat, path) 



