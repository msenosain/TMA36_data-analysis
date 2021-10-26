###############################################################################
# WES 
###############################################################################
## 1. Binary matrix
library(maftools)
library(dplyr)
library(tidyr)
maf_dir <- "data/TMA36_project/WES/processed/TwistWES_Tumor_pipeline.freq0.01.filtered.tsv_020921.maf"

mut = read.maf(maf = read.delim(maf_dir))
x = getClinicalData(mut)
x$Tumor_Sample_Barcode = as.character(x$Tumor_Sample_Barcode)
pt_ID <- sapply(strsplit(x$Tumor_Sample_Barcode, "pt"), "[[", 2)
pt_ID <- sapply(strsplit(pt_ID, "_"), "[[", 1)
CDE <- read.csv('data/TMA36_project/CDE/CDE_TMA36_2021SEPT21_DR_MF.csv')

CDE <- CDE[match(pt_ID, CDE$pt_ID),]
CDE <- cbind('Tumor_Sample_Barcode'=x$Tumor_Sample_Barcode, CDE)
mut = read.maf(maf = read.delim(maf_dir), clinicalData = CDE)

mut_df <- mut@data

h_sym <- unique(mut_df$Hugo_Symbol)
pt_ID <- unique(mut_df$Tumor_Sample_Barcode)

mut_dt <- matrix(0, nrow=length(h_sym), ncol=length(pt_ID))
rownames(mut_dt) <- h_sym
colnames(mut_dt) <- pt_ID


for (i in 1:length(h_sym)){
  k <- which(mut_df$Hugo_Symbol == h_sym[i])
  pt <- mut_df$Tumor_Sample_Barcode[k]
  mut_dt[i,pt] <- 1
}

colnames(mut_dt) <- sapply(strsplit(as.character(pt_ID), "pt"), "[[", 2)
colnames(mut_dt) <- sapply(strsplit(colnames(mut_dt), "_"), "[[", 1)
mut_dt <- t(mut_dt)

# Compute mutational load (number of mutations per patient)
#mut_load = data.frame(mut_load = apply(mut_dt, 1, sum))
mut_load <- data.frame(table(mut_df$Tumor_Sample_Barcode))
colnames(mut_load) <- c('pt_ID', 'mut_load')
mut_load$pt_ID <- sapply(strsplit(as.character(mut_load$pt_ID), "pt"), "[[", 2)
mut_load$pt_ID <- sapply(strsplit(as.character(mut_load$pt_ID), "_"), "[[", 1)


# Remove genes that are mutated in less than X% of samples
prcnt <- 0.1 # cutoff: genes must be mutated in >10 % of the samples
mut_dt_sum <- data.frame(t(mut_dt)) %>% 
  mutate(sum = rowSums(., na.rm = TRUE), 
         genes = rownames(.)) %>%
  arrange(., desc(sum)) %>%
  dplyr::select(genes, sum) %>%
  filter(., sum > nrow(mut_dt)*prcnt)

mut_dt <- data.frame(mut_dt[,mut_dt_sum$genes]) %>%
  mutate(pt_ID=rownames(.))%>%
  left_join(mut_load, ., by = 'pt_ID')

write.csv(mut_dt, file = "data/TMA36_project/WES/processed/wes_binary.csv")