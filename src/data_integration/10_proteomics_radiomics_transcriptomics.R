library(reticulate)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(cluster)
library(factoextra)
source('src/cytof/20_ClustAnnot_functions.R')

use_python("/Users/senosam/opt/anaconda3/envs/r-reticulate/bin/python", required = T)
py_config()
# Read data
cytof_freq <- read.csv("data/TMA36_project/CyTOF/processed/Data_paper2/both/cytof_freq.csv", row.names=1)
cytof_medianprot <- read.csv("data/TMA36_project/CyTOF/processed/Data_paper2/both/cytof_medianprot.csv", row.names=1)
rna_pathways <- read.csv("data/TMA36_project/RNA_Seq/processed/pathways_scores.csv", row.names=1)
rad_hm <- read.csv("data/TMA36_project/Radiomics/processed/rad_healthmyne.csv", row.names=1)
cde <- read.csv("data/TMA36_project/CDE/CDE_TMA36_2020FEB25_SA_MF.csv")

# get significant HM features with pairwise correlation
source_python('src/data_integration/pw_corr.py')
rad_hm_sigcor <- pw_corr()

# select significant variables for cytof and hm
cytof_freq_sig <- cytof_freq %>% select(., ECC_3,ECC_5,fmes_3,OtherI_4)
cytof_medianprot_sig <- cytof_medianprot %>% select(., HLA.DR)
rad_hm_sig <- rad_hm %>% select(rad_hm_sigcor$Y)

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
                            rho_cutoff=0.3,k_max=10){
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
  
  clusters
}

# Cluster features
set.seed(1)
#clusters_features <- cluster_by_corr(all_vars, adj_pval_cutoff=0.05, rho_cutoff=0.3, k_max=40)
k_feature = DetermineNumberOfClusters(t(all_vars_scaled), k_max = 15, ask_ft = F, arcsn_tr = F)
clusters_features <- data.frame('cluster'=kmeans(t(all_vars_scaled), k_feature, iter.max = 100)$cluster)

# Cluster patients
set.seed(1)
#clusters_patients <- cluster_by_corr(t(all_vars_scaled), adj_pval_cutoff=0.05, rho_cutoff=0, k_max=40)
k_patient = DetermineNumberOfClusters(all_vars_scaled, k_max = 15, ask_ft = F, arcsn_tr = F)
clusters_patients <- data.frame('cluster'=kmeans(all_vars_scaled, k_patient, iter.max = 100)$cluster)


# write cluster info into csv
write.csv(all_vars, file = 'data/TMA36_project/data_integration/cytof_rna_hm_raw.csv', row.names = TRUE)
write.csv(all_vars_scaled, file = 'data/TMA36_project/data_integration/cytof_rna_hm.csv', row.names = TRUE)
write.csv(clusters_patients, file = 'data/TMA36_project/data_integration/clusters_patients.csv', row.names = TRUE)
write.csv(clusters_features, file = 'data/TMA36_project/data_integration/clusters_features.csv', row.names = TRUE)

rad_hm_sigcor$`CI95%` <- as.character(rad_hm_sigcor$`CI95%`)
write.csv(rad_hm_sigcor, 'data/TMA36_project/Radiomics/processed/pwcorrelation_HMandSILA.csv', row.names = F)