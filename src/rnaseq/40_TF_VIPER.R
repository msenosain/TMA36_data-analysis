library(viper)
library(dorothea)
library(dplyr)
library(tidyverse)

msviper_mrs <- function(ls_preprocessed, data_collapsed, test_name, 
                        ref_name, file_name){
  
  # Dorothea to viper function
  dorothea2viper_regulons <- function(df) {
    regulon_list <- split(df, df$tf)
    viper_regulons <- lapply(regulon_list, function(regulon) {
      tfmode <- stats::setNames(regulon$mor, regulon$target)
      list(tfmode = tfmode, likelihood = rep(1, length(tfmode)))
    })
    return(viper_regulons)
  }
  
  # Generating regulons data
  data("dorothea_hs", package = "dorothea")
  regulons <- dorothea_hs %>%
    filter(confidence %in% c("A", "B"))
  regu <- dorothea2viper_regulons(regulons)
  
  # Generating test and ref data
  test_i <- which(ls_preprocessed$CDE$n_op2 %in% test_name)
  ref_i <- which(ls_preprocessed$CDE$n_op2 %in% ref_name)
  
  mat_test <- as.matrix(data_collapsed[,test_i])
  mat_ref <- as.matrix(data_collapsed[,ref_i])
  
  # Generating NULL model (test, reference)
  dnull <- ttestNull(mat_test, mat_ref, per=1000)
  
  # Generating signature
  signature <- rowTtest(mat_test, mat_ref)
  signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) * sign(signature$statistic))[, 1]
  
  # Running msVIPER
  mrs <- msviper(signature, regu, dnull, verbose = FALSE)
  mrs_table <- summary(mrs)
  
  
  regulons_mrs <- left_join(mrs_table, regulons, by = c('Regulon' = 'tf'))
  
  write.csv(mrs_table, file.path('data/TMA36_project/RNA_Seq/processed', paste0(file_name, '_MRS.csv')),
              row.names = FALSE, quote = FALSE)
  
  write.table(regulons_mrs, file.path('data/TMA36_project/RNA_Seq/processed', paste0(file_name, '.txt')) , 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
}

load('data/TMA36_project/RNA_Seq/processed/DE_analysis.RData')
result <- WGCNA::collapseRows(ls_preprocessed$vsd_mat,
                              rowGroup=ls_preprocessed$raw_counts$Feature_gene_name,
                              rowID=rownames(ls_preprocessed$vsd_mat),
                              method="MaxMean")

data <- data.frame(result$datETcollapsed)

# Indolent vs Aggressive  
msviper_mrs(ls_preprocessed, 
            data_collapsed = data, 
            test_name = 'agg', 
            ref_name = 'ind', 
            file_name = 'TF_VIPER_indvsagg')

# Indolent vs Intermediate  
msviper_mrs(ls_preprocessed, 
            data_collapsed = data, 
            test_name = 'int', 
            ref_name = 'ind', 
            file_name = 'TF_VIPER_indvsint')

# Intermediate vs Aggressive
msviper_mrs(ls_preprocessed, 
            data_collapsed = data, 
            test_name = 'agg', 
            ref_name = 'int', 
            file_name = 'TF_VIPER_intvsagg')

# Indolent vs Intermediate + Aggressive  
msviper_mrs(ls_preprocessed, 
            data_collapsed = data, 
            test_name = c('int','agg'),
            ref_name = 'ind', 
            file_name = 'TF_VIPER_indvsintagg')

# Intermediate vs Indolent + Aggressive
msviper_mrs(ls_preprocessed, 
            data_collapsed = data, 
            test_name = c('ind','agg'), 
            ref_name = 'int', 
            file_name = 'TF_VIPER_intvsindagg')

# Aggressive vs Indolent + Intermediate
msviper_mrs(ls_preprocessed, 
            data_collapsed = data, 
            test_name = c('ind','int'), 
            ref_name = 'agg', 
            file_name = 'TF_VIPER_aggvsindint')


# Indolent + Intermediate vs Aggressive
msviper_mrs(ls_preprocessed, 
            data_collapsed = data, 
            test_name = 'agg', 
            ref_name = c('ind','int'), 
            file_name = 'TF_VIPER_indintvsagg')


#############################################################
# TF by sample
#############################################################
# Generating regulons data
data("dorothea_hs", package = "dorothea")
regulons <- dorothea_hs %>%
  filter(confidence %in% c("A", "B"))
regu <- dorothea2viper_regulons(regulons)

vpres_4 <- t(viper(data, regu, verbose = FALSE, minsize = 4))
vpres_4 <- cbind('pt_ID' = ls_preprocessed$batch_info$pt_ID, vpres_4)

write.csv(vpres_4, file.path('data/TMA36_project/RNA_Seq/processed', 'TF_VIPER_all.csv'),
          row.names = FALSE, quote = FALSE)

