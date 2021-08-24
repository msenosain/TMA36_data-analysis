source('src/rnaseq/30_DEGanalysis.R')
source('src/rnaseq/30_DEGanalysis.R')
source('src/cytof/31_supervised_analysis_viz.R')
source('src/cytof/20_ClustAnnot_functions.R')
environment_set()

# Read in data
raw_counts <- read.csv('data/TMA36_project/RNA_Seq/processed/rawcounts_hg19.csv')
batch_info <- read.csv('data/TMA36_project/RNA_Seq/processed/rnaseq_batchinfo.csv')
CDE <- read.csv(file = 'data/TMA36_project/CDE/CDE_TMA36_2021JAN13_SA_MF.csv')

# Preprocessing data
ls_preprocessed <- preprocess_rna(raw_counts, batch_info, CDE,
                                  correct_batch = F, lowvargenesrm = T, prot_coding_only = F, xychr_rm = T)

# Add new columns for new labels
ind_idx <- which(ls_preprocessed$CDE$n_op2 == 'ind')
int_idx <- which(ls_preprocessed$CDE$n_op2 == 'int')
agg_idx <- which(ls_preprocessed$CDE$n_op2 == 'agg')

## ind vs int + agg
ls_preprocessed$CDE['indvs_int_agg'] <- ls_preprocessed$CDE$n_op2
ls_preprocessed$CDE$indvs_int_agg[c(int_idx, agg_idx)] <- 'int_agg'

## int vs ind + agg
ls_preprocessed$CDE['intvs_ind_agg'] <- ls_preprocessed$CDE$n_op2
ls_preprocessed$CDE$intvs_ind_agg[c(ind_idx, agg_idx)] <- 'ind_agg'

## agg vs ind + int
ls_preprocessed$CDE['aggvs_ind_int'] <- ls_preprocessed$CDE$n_op2
ls_preprocessed$CDE$aggvs_ind_int[c(ind_idx, int_idx)] <- 'ind_int'

#-------- INDOLENT VS AGGRESSIVE --------#
# DEG analysis
DE_res_indvsagg <- DE_analysis(ls_preprocessed, 
                      GeneBased=FALSE, 
                      pDataBased=TRUE,
                      NewCondition=FALSE,
                      cond_nm='n_op2',
                      two_levels=c('agg','ind'),
                      reference = 'ind')
# GSEA analysis
fgsea_res_indvsagg <- fgsea_analysis(DE_res_indvsagg)


#-------- INDOLENT VS INTERMEDIATE --------#
# DEG analysis
DE_res_indvsint <- DE_analysis(ls_preprocessed, 
                               GeneBased=FALSE, 
                               pDataBased=TRUE,
                               NewCondition=FALSE,
                               cond_nm='n_op2',
                               two_levels=c('int','ind'),
                               reference = 'ind')
# GSEA analysis
fgsea_res_indvsint <- fgsea_analysis(DE_res_indvsint)


#-------- INTERMEDIATE VS AGGRESSIVE --------#
# DEG analysis
DE_res_intvsagg <- DE_analysis(ls_preprocessed, 
                               GeneBased=FALSE, 
                               pDataBased=TRUE,
                               NewCondition=FALSE,
                               cond_nm='n_op2',
                               two_levels=c('agg','int'),
                               reference = 'int')
# GSEA analysis
fgsea_res_intvsagg <- fgsea_analysis(DE_res_intvsagg)


#-------- INDOLENT VS INT+AGGRESSIVE --------#
# DEG analysis
DE_res_indvsintagg <- DE_analysis(ls_preprocessed, 
                               GeneBased=FALSE, 
                               pDataBased=TRUE,
                               NewCondition=FALSE,
                               cond_nm='indvs_int_agg',
                               two_levels=c('int_agg','ind'),
                               reference = 'ind')
# GSEA analysis
fgsea_res_indvsintagg <- fgsea_analysis(DE_res_indvsintagg)


#-------- INT VS INDOLENT+AGGRESSIVE --------#
# DEG analysis
DE_res_intvsindagg <- DE_analysis(ls_preprocessed, 
                                  GeneBased=FALSE, 
                                  pDataBased=TRUE,
                                  NewCondition=FALSE,
                                  cond_nm='intvs_ind_agg',
                                  two_levels=c('ind_agg','int'),
                                  reference = 'int')
# GSEA analysis
fgsea_res_intvsindagg <- fgsea_analysis(DE_res_intvsindagg)


#-------- AGGRESSIVE VS INDOLENT+INT --------#
# DEG analysis
DE_res_aggvsindint <- DE_analysis(ls_preprocessed, 
                                  GeneBased=FALSE, 
                                  pDataBased=TRUE,
                                  NewCondition=FALSE,
                                  cond_nm='aggvs_ind_int',
                                  two_levels=c('ind_int','agg'),
                                  reference = 'agg')
# GSEA analysis
fgsea_res_aggvsindint <- fgsea_analysis(DE_res_aggvsindint)



###################################
# SAVE DATA
###################################
save(ls_preprocessed, 
     DE_res_indvsagg, fgsea_res_indvsagg,
     DE_res_indvsint, fgsea_res_indvsint,
     DE_res_intvsagg, fgsea_res_intvsagg,
     DE_res_indvsintagg, fgsea_res_indvsintagg,
     DE_res_intvsindagg, fgsea_res_intvsindagg,
     DE_res_aggvsindint, fgsea_res_aggvsindint, file = 'data/TMA36_project/RNA_Seq/processed/DE_analysis.RData')

