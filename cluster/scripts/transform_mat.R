taskID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

library(pacman)
p_load(ChIPseeker, TxDb.Hsapiens.UCSC.hg19.knownGene, BSgenome.Hsapiens.UCSC.hg19, clusterProfiler, rtracklayer,
GenomicRanges, ENCODExplorer, ggplot2, dplyr, reshape, heatmaps, abind, preprocessCore)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
genome = BSgenome.Hsapiens.UCSC.hg19

histone_marks = c("H3K4me3", "H3K36me3", "H3K4me1", "H3K27me3", "H3K9ac", "H2AFZ", "H4K20me1", "H3K9me3", "H3K4me2", "H3K27ac", "H3K79me2", "H3K9me1")
cell_lines = c("K562", "HepG2", "A549", "GM12878")
ENCODE_promoters_ext = readRDS("../data/ENCODE_promoters_ext.rds")
##############################
cl = readRDS("../data/hmm_training/accession_df.rds")
hm = readRDS("../data/promoters_HM.rds")
#ENCODE_promoters_ext = ENCODE_promoters_ext[taskID]
cl_subset = subset(cl, cell_type==cell_lines[taskID])

#### read multiple RDS files into a list 
pattern = paste0(cell_lines[taskID], "*")
files = list.files(path = "../data/hmm_training/matrices/", pattern = pattern)
files = paste("../data/hmm_training/matrices/", files, sep="")
files_id = as.numeric(unlist(lapply(strsplit(files, "[^0-9]+"), function(x) x=x[length(x)])))
names(files) = files_id
files_id = as.character(sort(files_id))
files_sorted = c()
for(i in 1:length(files_id)){files_sorted[i]=files[files_id[i]]}

binned_list = list()
for(i in 1:length(files_sorted)){binned_list[[i]] = readRDS(files_sorted[i])}
names(binned_list) = cl_subset$accession


# list of matrices to 3D array 
binned_list = abind(binned_list, along = 3) 
dimnames(binned_list)[[3]] = cl_subset$target
binned_data = lapply(seq_along(binned_list[,1,1]), function(x) return(binned_list[x,,]))

names = paste0(names(ENCODE_promoters_ext[taskID]), "_", c(1:length(ENCODE_promoters_ext[[taskID]])))
names(binned_data) = names

# reduce the training dataset to only the 'best' sample 
cl_spl = split(cl_subset, cl_subset$target)
binned_data = lapply(binned_data, function(z){colnames(z) = cl_subset$accession 
              return(z)})

.remove_outliers = function(x){
   names = c()
   for(i in 1:length(x)){
      #sample 100 seqs 
      spl = sample(binned_data, 500)
      spl = lapply(spl, function(y) y=y[, x[[i]]$accession])
      if(dim(x[[i]])[1]<2){
         names[i] = x[[i]]$accession
         next
      } else {
         v = lapply(spl, function(x) x=abs(var(x)))
         v = apply(sapply(v, function(x) apply(x, 1, mean)), 1, mean)
         names[i] = names(which.min(v)) 
      }
   }
   colnames=colnames(binned_data[[1]])[colnames(binned_data[[1]]) %in% names]
   binned_data = lapply(binned_data, function(x) x=x[, colnames])   
   binned_data = lapply(binned_data, function(x){colnames(x)=cl_subset$target[cl_subset$accession %in% colnames]
                                   return(x)})
return(binned_data)
}

binned_data = .remove_outliers(cl_spl)

# remove duplicates from training data 
ENCODE_promoters = readRDS("../data/ENCODE_promoters.rds")
ENCODE_promoters = ENCODE_promoters[-5]
names(ENCODE_promoters) = cell_lines
ENCODE_promoters = ENCODE_promoters[[cell_lines[taskID]]]

dupl = which(duplicated(ENCODE_promoters$geneId)) 
print(head(dupl), 20)
binned_data = binned_data[-dupl]
print("::::::")
print(binned_data[[1]][1:5, 1:10])

# normalisation 
.qnorm = function(x){
   colsum = apply(x, 2, sum)
   return(colsum)
}
colsum = t(sapply(binned_data, .qnorm))
colsum_norm = normalize.quantiles(colsum)
norm_factor = colsum/colsum_norm

data_norm = lapply(seq_along(binned_data), function(x) round(binned_data[[x]]/norm_factor[x,], 1))
# remove matrices with NaN values
x=sapply(data_norm, function(x) any(is.nan(x)))
data_norm = data_norm[!x]


# reduce the number of bins (1kb upstream/downstream)
binned_data = lapply(binned_data, function(x) x=x[c(16:35),])

name = names(ENCODE_promoters_ext[taskID])
path = paste0("../data/hmm_training/binned_hm_data_", name, ".rds") 
saveRDS(binned_data, path)





