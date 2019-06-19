#!/usr/bin/Rscript
#SBATCH -c 4
#SBATCH --mem 32G

cell_line = as.character(commandArgs(trailingOnly=TRUE))
library(pacman)
p_load(STAN, GenomicRanges, gplots, cowplot)


### import data
methods = list("NegativeBinomial", "PoissonLogNormal", "Bernoulli")
path = paste0("../data/hmm_training/binned_hm_data_", cell_line, ".rds")
training_data = readRDS(path)

### calculate average signal + produce heatmaps
.getHeatmap = function(x){
path = paste0("../data/HMMs/", x, "/", cell_line, "_viterbi.rds")
model_Viterbi = readRDS(path)
avg_signal = getAvgSignal(model_Viterbi, training_data)
 
## specify color palette
library(gplots)
heat = c("dark blue", "dodgerblue4", "darkred", "red", "orange", "gold", "yellow")
colfct = colorRampPalette(heat)
colpal_statemeans = colfct(25)

## define state order and colors
ord_nb = order(apply(avg_signal,1,max), decreasing=TRUE)
nStates = 24
statecols_nb = rainbow(nStates)
names(statecols_nb) = ord_nb

hm = heatmap.2(log(avg_signal+1)[as.character(ord_nb),], margins=c(8,7), srtCol=45,
        RowSideColors=statecols_nb[as.character(ord_nb)], dendrogram="none",
        Rowv=FALSE, Colv=FALSE, col=colpal_statemeans, trace="none",
        cellnote=round(avg_signal,1)[as.character(ord_nb),], notecol="black")

path = paste0("../data/HMMs/heatmaps/", cell_line, "_Bernoulli.png")
png(path)
dev.off()
}
#hm_list = lapply(methods, .getHeatmap)
.getHeatmap("Bernoulli")



### plot a list of heatmaps
#hms = plot_grid(hm_list, ncol=3, labels=methods)
#path = paste0("../data/HMMs/heatmaps/", cell_line, ".png")
#saveRDS(hms, path)


