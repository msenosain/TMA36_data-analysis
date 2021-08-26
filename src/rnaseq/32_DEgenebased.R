#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Gene symbol must be provided", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = paste0('data/TMA36_project/RNA_Seq/extra/', args[1], '_DE.Rdata')
}

# Load data
load('data/TMA36_project/RNA_Seq/processed/DE_lspreprocessed.RData')

gene_symbol <- args[1]
gn <- as.character(ls_preprocessed$raw_counts$Feature[which(ls_preprocessed$raw_counts$Feature_gene_name ==gene_symbol)])
if(length(gn)==0){
  stop(paste("Gene symbol", args[1], "does not exist!"), call.=FALSE)
}

source('src/rnaseq/30_DEGanalysis.R') 
environment_set()

DE_res <- DE_analysis(ls_preprocessed, 
                      GeneBased=TRUE, 
                      pDataBased=FALSE,
                      NewCondition=FALSE,
                      cond_nm= gn,
                      reference = 'low', 
                      #correct_gender=TRUE,
                      extremes_only=TRUE)

fgsea_res <- fgsea_analysis(DE_res)

save(gene_symbol, gn, DE_res, fgsea_res, file = args[2])