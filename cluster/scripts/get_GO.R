#!/usr/bin/Rscript 
#SBATCH -c 4
#SBATCH --mem 32G
#SBATCH --array 1-10

taskID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
arg = commandArgs(trailingOnly=TRUE)
cell_line = as.character(arg)
library(pacman)
p_load(STAN, GenomicRanges, ggplot2, reshape, dplyr, cowplot, gplots, 
       TxDb.Hsapiens.UCSC.hg19.knownGene, clusterProfiler, org.Hs.eg.db, preprocessCore, annotate)
ENCODE_promoters = readRDS("../data/ENCODE_promoters.rds")
cell_lines = c("K562", "HepG2", "A549", "GM12878")
names(ENCODE_promoters) = cell_lines
optimal_states = list(c(15, 10), c(15, 10), c(15, 10), c(15, 10))
names(optimal_states) = cell_lines

#######################
### GO ANNOTATION ###
#######################

# for each promoter cluster, determine and compare GO terms
.analyseGO = function(cell_line){
# load data and produce matrix   
path = paste0("../data/HMMs/Bernoulli/", cell_line, "_", optimal_states[[cell_line]][1], "_viterbi.rds")
viterbi_hmm = readRDS(path)   
viterbi_mat = do.call(rbind, viterbi_hmm)
clust = kmeans(viterbi_mat, centers = optimal_states[[cell_line]][2], iter.max = 15)

#partition seqs into groups
seq_id = row.names(viterbi_mat)
groups = unlist(lapply(seq_id, function(x) clust$cluster[[x]]))
group_df = data.frame(seq_id, groups)
groups = unique(group_df$group)

# GO profile for each group of promoters (arg = group)
# txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
all_genes = unique(ENCODE_promoters[[cell_line]]$geneId)

for(i in groups[taskID]){
   seqs_by_group = as.integer((row.names(subset(group_df, groups==i))))
   genes = unique(ENCODE_promoters[[cell_line]][seqs_by_group]$geneId)
   
   # cellular component
   ego1 = enrichGO(gene=genes, OrgDb = org.Hs.eg.db, universe = all_genes,  ont = "CC", keyType = 'ENTREZID', 
                   pAdjustMethod = "BH", pvalueCutoff  = 0.09, qvalueCutoff  = 0.05, readable = TRUE)
   
   path = paste0("../data/HMMs/GO_analysis/", cell_line,  "/CC_group_", i, ".png") 
   png(path, height=20, width=40, units="cm", res=200)
   bp = barplot(ego1, showCategory=12)
   dp = dotplot(ego1)
   print(plot_grid(bp, dp, ncol=2, labels = "AUTO"))
   dev.off()
   print("CC DONE")
   
   # molecular function 
   ego2 = enrichGO(gene=genes, OrgDb = org.Hs.eg.db, universe = all_genes,  
                   ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE)
   path = paste0("../data/HMMs/GO_analysis/", cell_line,  "/MF_group_", i, ".png") 
   png(path, height=20, width=35, units="cm", res=200)
   bp = barplot(ego2, showCategory=10)
   dp = dotplot(ego2)
   print(plot_grid(bp, dp, ncol=2, labels = "AUTO"))
   dev.off()
   print("MF DONE")
   
   # biological process
   ego3 = enrichGO(gene=genes, OrgDb = org.Hs.eg.db, universe = all_genes,  ont = "BP", 
                   pAdjustMethod = "BH", pvalueCutoff  = 0.1, qvalueCutoff  = 0.05, readable = TRUE)
   path = paste0("../data/HMMs/GO_analysis/", cell_line,  "/BP_group_", i, ".png") 
   png(path, height=20, width=25, units="cm", res=200)
   bp = barplot(ego3, showCategory=10)
   dp = dotplot(ego3)
   print(plot_grid(bp, dp, nrow=2, labels="AUTO"))
   dev.off()
   print("BP done")
       }
print(paste0("cell line analysed: ", cell_line))
}

.analyseGO(cell_line)


