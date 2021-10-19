###################################################################
# PATHWAYS SCORES USING GSVA
###################################################################
#https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_03_gsva.html

library(msigdbr)
library(GSVAdata)
library(GSVA)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(tidyverse)


# Load DE analysis data
load("data/TMA36_project/RNA_Seq/processed/DE_analysis.RData")

# Normalized gene expression data
vst_df <- ls_preprocessed$vsd_mat %>%
  as.data.frame(.)%>%
  mutate(ensembl_id=sapply(strsplit(rownames(.), "\\."), "[[", 1)) # remove variants info

# Collapse rows with simplified ENSEMBL ids
result <- WGCNA::collapseRows(vst_df[,-ncol(vst_df)],
                              rowGroup=vst_df$ensembl_id,
                              rowID=rownames(vst_df),
                              method="MaxMean")

vst_df <- data.frame(result$datETcollapsed) %>%
  mutate(ensembl_id=rownames(.))

# Map to ENTREZ ids
mapped_df <- data.frame("entrez_id"=mapIds(org.Hs.eg.db, 
                                           keys = vst_df$ensembl_id,
                                           keytype = "ENSEMBL",
                                           column = "ENTREZID",
                                           multiVals = "first")
                        ) %>%
  dplyr::filter(!is.na(entrez_id)) %>%
  tibble::rownames_to_column("Ensembl") %>%
  dplyr::inner_join(vst_df, by = c("Ensembl" = "ensembl_id")) %>%
  tibble::column_to_rownames('Ensembl')

# Remove duplicated ENTREZ ids
result <- WGCNA::collapseRows(mapped_df[,-1],
                              rowGroup=mapped_df$entrez_id,
                              rowID=rownames(mapped_df),
                              method="MaxMean")

filtered_mapped_df <- as.matrix(result$datETcollapsed)



# Get fgsea results for Hallmark
pathways_hallmark <- rbind(fgsea_res_indvsagg$res_hm, 
                     fgsea_res_indvsint$res_hm, 
                     fgsea_res_intvsagg$res_hm)%>%
  filter(., padj<0.05,abs(NES) >=1.5) %>% # Select only significant pathways
  distinct(pathway) %>%
  pull()

hm_geneset <- msigdbr::msigdbr(
  species = "Homo sapiens", 
  category = "H") %>%
  filter(gs_name %in% pathways_hallmark)

hallmark_list <- split(
  hm_geneset$entrez_gene, # The genes we want split into pathways
  hm_geneset$gs_name # The pathways made as the higher levels of the list
)

# Get fgsea results for REACTOME
pathways_reactome <- rbind(fgsea_res_indvsagg$res_rtm,
                           fgsea_res_indvsint$res_rtm,
                           fgsea_res_intvsagg$res_rtm)%>%
  filter(., padj<0.05,abs(NES) >=1.5) %>%
  distinct(pathway, .keep_all = T)

reactome_list <- str_split(pathways_reactome$leadingEdge, pattern = ' ')
names(reactome_list) <- pathways_reactome$pathway

gene_set_list <- c(hallmark_list, reactome_list)

# Perform GSVA
gsva_results <- gsva(
  filtered_mapped_df,
  gene_set_list,
  method = "gsva",
  # Appropriate for our vst transformed data
  kcdf = "Gaussian",
  # Minimum gene set size
  min.sz = 15, #1
  # Maximum gene set size
  max.sz = 500,
  # Compute Gaussian-distributed scores
  mx.diff = TRUE,
  # Don't print out the progress bar
  verbose = FALSE
)

colnames(gsva_results) <- ls_preprocessed$batch_info$pt_ID
gsva_results <- t(gsva_results)
write.csv(gsva_results, file = 'data/TMA36_project/RNA_Seq/processed/pathways_gsva.csv')

