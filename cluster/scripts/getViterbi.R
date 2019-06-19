#!/usr/bin/Rscript
#SBATCH -c 4 
#SBATCH --mem 32G
#SBATCH --array 1-5

taskID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
cell_line = as.character(commandArgs(trailingOnly=TRUE))

### import packages and data
library(pacman)
p_load(STAN, GenomicRanges)
ENCODE_promoters_ext = readRDS("../data/ENCODE_promoters_ext.rds")

### get models corresponding to cell line + parallelise by method  
methods = c("Bernoulli", "PoissonLogNormal", "NegativeBinomial")
nStates = c(10, 15, 20, 25, 30)
pattern = paste0(cell_line, "_", methods[1], "_", nStates[taskID], "_HMM.rds") 
model_name = list.files(path = "../data/HMMs/tmpdir/", pattern = pattern)
path = paste0("../data/HMMs/tmpdir/", model_name)
model = readRDS(path)

### get training dataset
path = paste0("../data/hmm_training/binned_hm_data_", cell_line, ".rds")
training_data = readRDS(path)

# binarise if Bernoulli and apply Viterbi algorithm to annotate the bins
if(methods[1]=="Bernoulli"){
training_data = binarizeData(training_data, thresh = 1e-04)
model_Viterbi = getViterbi(model, training_data)
} else {
cell_type_lst = list(cell_type=grep(cell_line, names(training_data)))
size_factors = getSizeFactors(training_data, cell_type_lst)
model_Viterbi = getViterbi(model, training_data, size_factors)
}
path = paste0("../data/HMMs/", methods[1], "/", cell_line, "_", nStates[taskID], "_viterbi.rds")
saveRDS(model_Viterbi, path)

### convert to GRanges
model_ranges = viterbi2GRanges(model_Viterbi, ENCODE_promoters_ext[[cell_line]], binSize=100)
path = paste0("../data/HMMs/", methods[1], "/", cell_line,"_", nStates[taskID], "_ranges", ".rds")
saveRDS(model_ranges, path)

