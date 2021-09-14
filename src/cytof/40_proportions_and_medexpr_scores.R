###############################################################################
# CyTOF
###############################################################################

source("src/cytof/31_supervised_analysis_viz.R")
source("src/cytof/20_ClustAnnot_functions.R")

# Load data
load("data/TMA36_project/CyTOF/processed/Data_paper2/both/cellclusters.RData")
ref <- ref[-42,] # sample 13376 has been collected twice and will be merged, only one ref is needed
load("data/TMA36_project/CyTOF/processed/Data_paper2/both/percent_pt.RData")
load("data/TMA36_project/CyTOF/processed/Data_paper2/both/protein_correlations.RData")

#--------------------- Get Frequencies table ---------------------#

# Epithelial, Stroma, Immune, Endothelial, Mesenchymal, Immune subtypes, clusters
freq_cytof <- cbind('Epithelial'= prcnt_celltypes$Epithelial,
                    'Stroma' = prcnt_celltypes$Fib_Mesenchymal+prcnt_celltypes$Endothelial,
                    'Endothelial' = prcnt_celltypes$Endothelial,
                    'Fib_Mesenchymal' = prcnt_celltypes$Fib_Mesenchymal,
                    'Immune' = prcnt_celltypes$Immune,
                    prcnt_subtypes[,4:8],
                    prcnt_clusters)

# center and scale
freq_cytof <- scale(freq_cytof)

# write
write.csv(freq_cytof, file = 'data/TMA36_project/CyTOF/processed/Data_paper2/both/cytof_freq.csv')

#--------------------- Get Median bulk protein expression table ---------------------#

ft_cols <- c(15,17:31, 33:35, 37:51)
df_cluster <- denoisingCTF::t_asinh(cbind(annot_df[,ft_cols], 
                                          'pt_ID'= annot_df[,'pt_ID']))
cl_median <- aggregate(df_cluster[,-ncol(df_cluster)], list(df_cluster$pt_ID), median)
rownames(cl_median)<- cl_median$Group.1
cl_median <- cl_median[,-1]
colnames(cl_median) <- gsub(".*_",  "",colnames(cl_median)) 

# center and scale
cl_median <- scale(cl_median)

# write
write.csv(cl_median, file = 'data/TMA36_project/CyTOF/processed/Data_paper2/both/cytof_medianprot.csv')
