library(reticulate)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(cluster)
library(factoextra)
library(RCy3)
library(igraph)
source('src/cytof/20_ClustAnnot_functions.R')

use_python("/Users/senosam/opt/anaconda3/envs/r-reticulate/bin/python", required = T)
py_config()
# Read data
cytof_freq <- read.csv("data/TMA36_project/CyTOF/processed/Data_paper2/both/cytof_freq.csv", row.names=1)
cytof_medianprot <- read.csv("data/TMA36_project/CyTOF/processed/Data_paper2/both/cytof_medianprot.csv", row.names=1)
#rna_pathways <- read.csv("data/TMA36_project/RNA_Seq/processed/pathways_scores.csv", row.names=1)
rna_pathways <- read.csv("data/TMA36_project/RNA_Seq/processed/pathways_gsva15.csv", row.names=1)
rad_hm <- read.csv("data/TMA36_project/Radiomics/processed/rad_healthmyne.csv", row.names=1)
cde <- read.csv("data/TMA36_project/CDE/CDE_TMA36_2020FEB25_SA_MF.csv")

# get significant HM features with pairwise correlation
source_python('src/data_integration/pw_corr.py')
rad_hm_sigcor <- pw_corr()

# select significant variables for cytof and hm
cytof_freq_sig <- cytof_freq %>% dplyr::select(., ECC_3,ECC_5,fmes_3,OtherI_4)
cytof_medianprot_sig <- cytof_medianprot %>% dplyr::select(., HLA.DR)
rad_hm_sig <- rad_hm %>% dplyr::select(rad_hm_sigcor$Y)

# join significant vars from all datasets in one df
all_vars <- inner_join(rownames_to_column(cytof_freq_sig, var = 'pt_ID'), 
                       rownames_to_column(cytof_medianprot_sig, var = 'pt_ID'), 
                       by = 'pt_ID')
all_vars <- inner_join(all_vars, 
                       rownames_to_column(rad_hm_sig, var = 'pt_ID'), 
                       by = 'pt_ID')
all_vars <- inner_join(all_vars, 
                       rownames_to_column(rna_pathways, var = 'pt_ID'), 
                       by = 'pt_ID')
rownames(all_vars) <- all_vars$pt_ID
all_vars$pt_ID <- NULL
colnames(all_vars) = gsub('\\.', ' ', colnames(all_vars))

# scale and center
all_vars_scaled <- scale(all_vars)

# Match CDE
cde <- cde[match(rownames(all_vars_scaled), cde$pt_ID),]

# Add WES info to CDE
wes = read.csv("data/TMA36_project/WES/processed/wes_binary.csv")
colnames(wes)[1] = 'pt_ID'
cde = left_join(cde, wes, by = 'pt_ID')

#----------------------- Clustering-----------------------#

cluster_by_corr <- function(data, adj_pval_cutoff=0.05, 
                            rho_cutoff=0.3,k_max=10, corr_only = FALSE){
  # Compute similarity matrix
  res <- Hmisc::rcorr(as.matrix(data), type = 'spearman') 
  rho <- res$r
  
  corrected_pvals <- p.adjust(res$P, method = 'BH')
  corrected_pvals <- matrix(corrected_pvals, nrow = ncol(res$P), 
                            ncol = ncol(res$P))
  # Remove non-significant correlations
  rho[corrected_pvals > adj_pval_cutoff] <- 0
  
  # Remove low correlations
  rho[abs(rho)<rho_cutoff] <- 0
  
  if(corr_only){
    return(rho)
  }else{
    # Compute dissimilarity matrix
    dis_ft <- as.dist((1-rho)/2) # https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/speardis.htm
    
    # Finding the optimal k
    k_opt <- fviz_nbclust((1-rho)/2, pam, method = "wss", print.summary=T, k.max = k_max) #https://www.statology.org/k-medoids-in-r/
    k_opt <- findElbow(k_opt$data$y, k_max)
    
    # create k-medoids clustering with k clusters
    ft_kmedoids <- pam(dis_ft, k_opt) 
    clusters <- as.data.frame(ft_kmedoids$clustering)
    rownames(clusters) = gsub('\\.', ' ', rownames(clusters))
    colnames(clusters) = 'cluster'
    
    return(clusters)
  }
}

# K means clustering
set.seed(1)
#clusters_features <- cluster_by_corr(all_vars, adj_pval_cutoff=0.05, rho_cutoff=0.3, k_max=40)
k_feature = DetermineNumberOfClusters(t(all_vars_scaled), k_max = 20, ask_ft = F, arcsn_tr = F) # 15
clusters_features <- data.frame('cluster'=kmeans(t(all_vars_scaled), k_feature, iter.max = 100)$cluster)

set.seed(1)
#clusters_patients <- cluster_by_corr(t(all_vars_scaled), adj_pval_cutoff=0.05, rho_cutoff=0, k_max=40)
k_patient = DetermineNumberOfClusters(all_vars_scaled, k_max = 10, ask_ft = F, arcsn_tr = F) # 15
clusters_patients <- data.frame('cluster'=kmeans(all_vars_scaled, k_patient, iter.max = 100)$cluster)

# Hierarchical clustering 
k = findElbow(fviz_nbclust(t(all_vars_scaled), hcut, method = "wss", print.summary=T, k.max = 20)$data$y, 20)
clusters_features2 <- data.frame('cluster'=hcut(t(all_vars_scaled), k, iter.max = 100)$cluster)

k = findElbow(fviz_nbclust(all_vars_scaled, hcut, method = "wss", print.summary=T, k.max = 10)$data$y, 10)
clusters_patients2 <- data.frame('cluster'=hcut(all_vars_scaled, k, iter.max = 100)$cluster)

# write cluster info into csv
write.csv(all_vars, file = 'data/TMA36_project/data_integration/cytof_rna_hm_raw.csv', row.names = TRUE)
write.csv(all_vars_scaled, file = 'data/TMA36_project/data_integration/cytof_rna_hm.csv', row.names = TRUE)
write.csv(clusters_patients, file = 'data/TMA36_project/data_integration/clusters_patients.csv', row.names = TRUE)
write.csv(clusters_features, file = 'data/TMA36_project/data_integration/clusters_features.csv', row.names = TRUE)
write.csv(clusters_patients2, file = 'data/TMA36_project/data_integration/clusters_patients2.csv', row.names = TRUE)
write.csv(clusters_features2, file = 'data/TMA36_project/data_integration/clusters_features2.csv', row.names = TRUE)

rad_hm_sigcor$`CI95%` <- as.character(rad_hm_sigcor$`CI95%`)
write.csv(rad_hm_sigcor, 'data/TMA36_project/Radiomics/processed/pwcorrelation_HMandSILA.csv', row.names = F)


# NETWORKS

ft_mat = cluster_by_corr(all_vars_scaled, adj_pval_cutoff=0.05, rho_cutoff=0, corr_only = T)
ft_network <- graph_from_adjacency_matrix(ft_mat, weighted=T, mode="undirected", 
                                       diag=F)

pt_mat = cluster_by_corr(t(all_vars_scaled), adj_pval_cutoff=0.05, rho_cutoff=0, corr_only = T)
pt_network <- graph_from_adjacency_matrix(pt_mat, weighted=T, mode="undirected", 
                                          diag=F)

# Open Cytoscape and confirm connection
cytoscapePing()
createNetworkFromIgraph(ft_network,"ft_network")
createNetworkFromIgraph(pt_network,"pt_network")
loadTableData(clusters_features)
loadTableData(clusters_patients)
