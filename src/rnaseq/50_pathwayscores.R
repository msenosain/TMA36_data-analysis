library(dplyr)
library(tidyverse)
library(biomaRt)
library(org.Hs.eg.db)
library(annotate)

################################################################################
# Function to compute pathways scores
################################################################################
pathway_score <- function(pth_res_ls, dt, REACTOME=FALSE, scale_=TRUE){
  
  # filter pathways (padj < 0.05, NES >= 1.5)
  pth_res_all <- matrix(NA, nrow = 1, 
                        ncol = ncol(pth_res_ls[[1]]), 
                        dimnames = list(1,colnames(pth_res_ls[[1]])))
  
  for(p in 1:length(pth_res_ls)){
    pth_res_ls[[p]] <-  pth_res_ls[[p]] %>%
      filter(., padj < 0.05,
             abs(NES) >=1.5)
    pth_res_all <- rbind(pth_res_all, pth_res_ls[[p]])
    if(p==length(pth_res_ls)){
      pth_res_all <- pth_res_all[-1,]
      pth_res_all <- distinct(pth_res_all, pathway, .keep_all = TRUE)
      rownames(pth_res_all) <- NULL
    }
  }
  
  paths=pth_res_all$pathway
  
  # create a list with the names of the DE pathways and the leading edge genes
  pathexpl=list()
  for (i in 1:nrow(pth_res_all)) {
    gsym <- unlist(str_split(pth_res_all$leadingEdge[i], pattern = ' '))
    if(REACTOME){
      gsym <- getSYMBOL(gsym, data='org.Hs.eg')
    }
    pathexpl[[paths[i]]] <- colSums(dt[gsym,])
  }
  
  # Create a df form the list
  pathexpdf <- Reduce(rbind, pathexpl)
  rownames(pathexpdf) <- names(pathexpl)
  pathexpdf <- pathexpdf[complete.cases(pathexpdf),]
  
  # Center and scale
  if(scale_){
    pathexpdf <- data.frame(scale(t(pathexpdf)))
  }else{
    pathexpdf <- data.frame(t(pathexpdf))
  }
  colnames(pathexpdf) <- paths
  
  return(pathexpdf)
  
}


################################################################################
# Run function
################################################################################


# Load DE analysis data
load("data/TMA36_project/RNA_Seq/processed/DE_analysis.RData")

# Normalized gene expression data
result <- WGCNA::collapseRows(ls_preprocessed$vsd_mat,
                              rowGroup=ls_preprocessed$raw_counts$Feature_gene_name,
                              rowID=rownames(ls_preprocessed$vsd_mat),
                              method="MaxMean")

dt <- data.frame(result$datETcollapsed)

pthres_hm_ls <- list(fgsea_res_indvsagg$res_hm, 
                   fgsea_res_indvsint$res_hm, 
                   fgsea_res_intvsagg$res_hm)

pthres_rtm_ls <- list(fgsea_res_indvsagg$res_rtm, 
                     fgsea_res_indvsint$res_rtm, 
                     fgsea_res_intvsagg$res_rtm)

hm_pth_scores <- pathway_score(pthres_hm_ls, dt, REACTOME=FALSE)
rtm_pth_scores <- pathway_score(pthres_rtm_ls, dt, REACTOME=TRUE)

pth_scores <- cbind(hm_pth_scores, rtm_pth_scores)
row.names(pth_scores) <- ls_preprocessed$batch_info$pt_ID

# Write csv files
write.csv(pth_scores, file = 'data/TMA36_project/RNA_Seq/processed/pathways_scores.csv')


