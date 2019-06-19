#!/usr/bin/Rscript
#SBATCH -c 4
#SBATCH --mem 32G
#SBATCH --array 1-5
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
nStates = c(10, 15, 20, 25, 30)
maxIters = 500
method = "Bernoulli"
print(cell_line)
print(nStates[taskID])
sample_binary = binarizeData(obs=sample, thresh = 1e-04)
hmm =  initHMM(sample_binary, nStates[taskID], method)

# find optimal number of iterations
hmm_fitted = fitHMM(sample_binary, hmm, maxIters = maxIters)
path = paste0("../data/HMMs/tmpdir/", cell_line, "_", method, "_", nStates[taskID],"_HMM.rds")
saveRDS(hmm_fitted, path)

log_L = hmm_fitted@LogLik
iterations = c(1:length(log_L))
data = data.frame(iterations, log_L)

line_1 = ggplot(data = data, aes(x=iterations, y=log_L)) +
   geom_line()+
   geom_point(alpha=0.2, color = "red")+
   labs(x = "Number of iterations", y = "log Likelihood", title=paste0(cell_line, "_", method, "::", "all iterations"))
      
line_2 = ggplot(data = data, aes(x=iterations, y=log_L)) +
   geom_line()+
   geom_point(alpha=0.2, color = "red")+
   labs(x = "Number of iterations", y = "log Likelihood", title=paste0(cell_line, "_", method,"::","Iterations 1-10 removed"))+
   coord_cartesian(ylim = c(data$log_L[10], data$log_L[length(data$log_L)]+1000))
plot = plot_grid(line_1, line_2, ncol=2)

path = paste0("../data/HMMs/tmpdir/", cell_line, "_", method,"_", nStates[taskID], ".png")
ggsave(path)

