#!/usr/bin/Rscript
#SBATCH -c 4
#SBATCH --mem 32G

taskID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
arg = commandArgs(trailingOnly=TRUE)
cell_line = as.character(arg)

library(pacman)
p_load(STAN, GenomicRanges, ggplot2, reshape, dplyr, cowplot)
binned_K562 = readRDS("../data/hmm_training/binned_hm_data_K562.rds")
binned_HepG2 = readRDS("../data/hmm_training/binned_hm_data_HepG2.rds")
binned_A549 = readRDS("../data/hmm_training/binned_hm_data_A549.rds")
binned_GM12878 = readRDS("../data/hmm_training/binned_hm_data_GM12878.rds")
binned_data = list(binned_K562, binned_HepG2, binned_A549, binned_GM12878)
names(binned_data) = c("K562", "HepG2", "A549", "GM12878")
sample = binned_data[[cell_line]]  
nStates = 12
maxIters = 500
method = "NegativeBinomial"

celltype_lst = list(cell_type=grep(cell_line, names(binned_data[[cell_line]])))
size_factors = getSizeFactors(sample, celltype_lst)
hmm =  initHMM(sample, nStates, method, size_factors)

# find optimal number of iterations
hmm_fitted = fitHMM(sample, hmm, maxIters = maxIters)
path = paste0("../data/HMMs/tmpdir/", cell_line, "_", method, "_", "HMM.rds")
saveRDS(hmm_fitted, path)

log_L = hmm_fitted@LogLik
iterations = c(1:length(log_L))
data = data.frame(iterations, log_L)

line_1 = ggplot(data, aes(x=iterations, y=log_L)) +
   geom_line()+
   geom_point(alpha=0.2, color = "red")+
   labs(x = "Number of iterations", y = "log Likelihood", title= paste0(cell_line, "_", method, "::", "all iterations"))
   
line_2 = ggplot(data, aes(x=iterations, y=log_L)) +
   geom_line()+
   geom_point(alpha=0.2, color = "red")+
   labs(x = "Number of iterations", y = "log Likelihood", title=paste0(cell_line, "_", method, "::", "Iterations 1-10 removed")) +
   coord_cartesian(ylim = c(data$log_L[10], data$log_L[length(data$log_L)]+1000))
plot = plot_grid(line_1, line_2, ncol=2)

path = paste0("../data/HMMs/tmpdir/", cell_line, "_", method, ".png")
ggsave(path)



