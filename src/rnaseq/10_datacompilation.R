library(dplyr)

###############################################################################
# Hg19
###############################################################################

#-------------- First Batch --------------#

counts_3388YZ <- read.delim(file = "data/TMA36_project/RNA_Seq/raw/Shilin/Hg19/RnaSeq_3388YZ_Tumor/RnaSeq_3388YZ_Tumor.count.txt")
fpkm_3388YZ <- read.delim(file = "data/TMA36_project/RNA_Seq/raw/Shilin/Hg19/RnaSeq_3388YZ_Tumor/RnaSeq_3388YZ_Tumor.fpkm.txt")

p_3388YZ <- readxl::read_excel("data/TMA36_project/RNA_Seq/raw/Shilin/Vantage_info/TMA 36 RNA DNA.xlsx", 
                               sheet = "3388-YZ", col_types = c("text", 
                                                                "text"))
# Match ptID
x <- match(paste0('pt.', p_3388YZ$pt_ID), colnames(counts_3388YZ)[9:ncol(counts_3388YZ)]) + 8

counts_3388YZ <- counts_3388YZ[,c(1:8,x)]
fpkm_3388YZ <- fpkm_3388YZ[,c(1:8,x)]

table(colnames(counts_3388YZ)[9:70] == paste0('pt.', p_3388YZ$pt_ID)) # confirm
table(colnames(fpkm_3388YZ)[9:70] == paste0('pt.', p_3388YZ$pt_ID)) # confirm

colnames(counts_3388YZ)[9:70] = p_3388YZ$Vantage_ID
colnames(fpkm_3388YZ)[9:70] = p_3388YZ$Vantage_ID

p_3388YZ['Batch'] = 1

counts_3388YZ$Feature_gene_name1 <- NULL
fpkm_3388YZ$Feature_gene_name1 <- NULL

p_3388YZ[37, 2] <- '14301'
p_3388YZ[56, 2] <- '8356'

#-------------- Second Batch --------------#

counts_4163YZ <- read.delim(file = "data/TMA36_project/RNA_Seq/raw/Shilin/Hg19/RnaSeq_4163YZ_Tumor/RnaSeq_3388YZ_Tumor/RnaSeq_3388YZ_Tumor.count.txt")
fpkm_4163YZ <- read.delim(file = "data/TMA36_project/RNA_Seq/raw/Shilin/Hg19/RnaSeq_4163YZ_Tumor/RnaSeq_3388YZ_Tumor/RnaSeq_3388YZ_Tumor.fpkm.txt")

p_4163YZ <- readxl::read_excel("data/TMA36_project/RNA_Seq/raw/Shilin/Vantage_info/TMA 36 RNA DNA.xlsx", 
    sheet = "4163-YZ", col_types = c("text", 
        "text"))

x <- match(paste0('X',gsub('-', '.', p_4163YZ$Vantage_ID)), colnames(counts_4163YZ)[8:ncol(counts_4163YZ)]) + 7

counts_4163YZ <- counts_4163YZ[,c(1:7,x)]
fpkm_4163YZ <- fpkm_4163YZ[,c(1:7,x)]

table(colnames(counts_4163YZ)[8:36] == paste0('X',gsub('-', '.', p_4163YZ$Vantage_ID)))
table(colnames(fpkm_4163YZ)[8:36] == paste0('X',gsub('-', '.', p_4163YZ$Vantage_ID)))

colnames(counts_4163YZ)[8:36] = p_4163YZ$Vantage_ID
colnames(fpkm_4163YZ)[8:36] = p_4163YZ$Vantage_ID

p_4163YZ['Batch'] = 2

p_4163YZ <- p_4163YZ[-29,]
counts_4163YZ[,36] <- NULL
fpkm_4163YZ[,36] <- NULL

#-------------- Merge --------------#

counts_all = dplyr::inner_join(counts_3388YZ, counts_4163YZ, by = colnames(counts_3388YZ)[1:7])
fpkm_all = dplyr::inner_join(fpkm_3388YZ, fpkm_4163YZ, by = colnames(fpkm_3388YZ)[1:7])
p_all <- rbind(p_3388YZ, p_4163YZ)

# Edit names
p_all$Vantage_ID = paste0('R',gsub('-', '_', p_all$Vantage_ID))
colnames(counts_all)[8:ncol(counts_all)] = p_all$Vantage_ID
colnames(fpkm_all)[8:ncol(fpkm_all)] = p_all$Vantage_ID

# Normalizing by gene length (TPM)
tpm <- function(counts, lengths) {
  x <- counts / lengths
  tpm.mat <- t( t(x) * 1e6 / colSums(x) )
  tpm.mat
}

tpm_all <- cbind(counts_all[,1:7], tpm(counts_all[,8:ncol(counts_all)], counts_all$Feature_length))

#-------------- Remove low quality samples and duplicates --------------#

# Low quality samples
qc_R4163 <- paste0('R4163_YZ_', c('4','12','13','14','15','27'))
qc_R3388 <- paste0('R3388_YZ_', c('8', '45'))

# Remove duplicates
dupl <- c('11817', '12889', '12929', '15002') #11840 13034 12890 13155
# Vector of duplicates and lq samples
unique_vntg <- p_all %>%
  filter(., !(pt_ID %in% dupl &
                grepl("R4163",Vantage_ID)),
         !(Vantage_ID %in% c(qc_R4163,qc_R3388))) %>%
  dplyr::select(., Vantage_ID) %>%
  pull()

# Clean the merged data
keep_cols <- c(colnames(counts_all)[1:7], unique_vntg)
counts_all <- counts_all[keep_cols]
fpkm_all <- fpkm_all[keep_cols]
tpm_all <- tpm_all[keep_cols]
p_all <- p_all[which(p_all$Vantage_ID %in% unique_vntg),]

#-------------- Write data --------------#
write.csv(counts_all, "data/TMA36_project/RNA_Seq/processed/rawcounts_hg19.csv", row.names = FALSE)
write.csv(fpkm_all, "data/TMA36_project/RNA_Seq/processed/fpkm_hg19.csv", row.names = FALSE)
write.csv(tpm_all, "data/TMA36_project/RNA_Seq/processed/tpm_hg19.csv", row.names = FALSE)
write.csv(p_all, "data/TMA36_project/RNA_Seq/processed/rnaseq_batchinfo.csv", row.names = FALSE)

###############################################################################
# Hg38
###############################################################################

#-------------- First Batch --------------#

counts_3388YZ <- read.delim(file = "data/TMA36_project/RNA_Seq/raw/Shilin/Hg38/20190827_3388YZ_Tumor_RnaSeq/RnaSeq_3388YZ_Tumor_hg38/RnaSeq_3388YZ_Tumor_hg38.count")
fpkm_3388YZ <- read.delim(file = "data/TMA36_project/RNA_Seq/raw/Shilin/Hg38/20190827_3388YZ_Tumor_RnaSeq/RnaSeq_3388YZ_Tumor_hg38/RnaSeq_3388YZ_Tumor_hg38.fpkm.tsv")

p_3388YZ <- readxl::read_excel("data/TMA36_project/RNA_Seq/raw/Shilin/Vantage_info/TMA 36 RNA DNA.xlsx", 
                               sheet = "3388-YZ", col_types = c("text", 
                                                                "text"))
# Match ptID
x <- match(paste0('X',gsub('-', '.', p_3388YZ$Vantage_ID)), colnames(counts_3388YZ)[8:ncol(counts_3388YZ)]) + 7

counts_3388YZ <- counts_3388YZ[,c(1:7,x)]
fpkm_3388YZ <- fpkm_3388YZ[,c(1:7,x)]

colnames(counts_3388YZ)[8:ncol(counts_3388YZ)] = p_3388YZ$Vantage_ID
colnames(fpkm_3388YZ)[8:ncol(fpkm_3388YZ)] = p_3388YZ$Vantage_ID

p_3388YZ['Batch'] = 1
p_3388YZ[37, 2] <- '14301'
p_3388YZ[56, 2] <- '8356'

#-------------- Second Batch --------------#
counts_4163YZ <- read.delim(file = "data/TMA36_project/RNA_Seq/raw/Shilin/Hg38/20200211_4163YZ_TumorBatch2_RnaSeq/RnaSeq_4163YZ_Tumor_Hg38/RnaSeq_4163YZ_Tumor_Hg38.count")
fpkm_4163YZ <- read.delim(file = "data/TMA36_project/RNA_Seq/raw/Shilin/Hg38/20200211_4163YZ_TumorBatch2_RnaSeq/RnaSeq_4163YZ_Tumor_Hg38/RnaSeq_4163YZ_Tumor_Hg38.fpkm.tsv")

p_4163YZ <- readxl::read_excel("data/TMA36_project/RNA_Seq/raw/Shilin/Vantage_info/TMA 36 RNA DNA.xlsx", 
                               sheet = "4163-YZ", col_types = c("text", 
                                                                "text"))

# Match ptID
x <- match(paste0('X',gsub('-', '.', p_4163YZ$Vantage_ID)), colnames(counts_4163YZ)[8:ncol(counts_4163YZ)]) + 7

counts_4163YZ <- counts_4163YZ[,c(1:7,x)]
fpkm_4163YZ <- fpkm_4163YZ[,c(1:7,x)]

colnames(counts_4163YZ)[8:ncol(counts_4163YZ)] = p_4163YZ$Vantage_ID
colnames(fpkm_4163YZ)[8:ncol(fpkm_4163YZ)] = p_4163YZ$Vantage_ID

p_4163YZ['Batch'] = 2

p_4163YZ <- p_4163YZ[-29,]
counts_4163YZ[,36] <- NULL
fpkm_4163YZ[,36] <- NULL

#-------------- Merge --------------#

counts_all = dplyr::inner_join(counts_3388YZ, counts_4163YZ, by = colnames(counts_3388YZ)[1:7])
fpkm_all = dplyr::inner_join(fpkm_3388YZ, fpkm_4163YZ, by = colnames(fpkm_3388YZ)[1:7])
p_all <- rbind(p_3388YZ, p_4163YZ)

# Edit names
p_all$Vantage_ID = paste0('R',gsub('-', '_', p_all$Vantage_ID))
colnames(counts_all)[8:ncol(counts_all)] = p_all$Vantage_ID
colnames(fpkm_all)[8:ncol(fpkm_all)] = p_all$Vantage_ID

# Normalizing by gene length (TPM)
tpm <- function(counts, lengths) {
  x <- counts / lengths
  tpm.mat <- t( t(x) * 1e6 / colSums(x) )
  tpm.mat
}

tpm_all <- cbind(counts_all[,1:7], tpm(counts_all[,8:ncol(counts_all)], counts_all$Feature_length))

#-------------- Remove low quality samples and duplicates --------------#

# Low quality samples
qc_R4163 <- paste0('R4163_YZ_', c('4','12','13','14','15','27'))
qc_R3388 <- paste0('R3388_YZ_', c('8', '45'))

# Remove duplicates
dupl <- c('11817', '12889', '12929', '15002') #11840 13034 12890 13155
# Vector of duplicates and lq samples
unique_vntg <- p_all %>%
  filter(., !(pt_ID %in% dupl &
                grepl("R4163",Vantage_ID)),
         !(Vantage_ID %in% c(qc_R4163,qc_R3388))) %>%
  dplyr::select(., Vantage_ID) %>%
  pull()

# Clean the merged data
keep_cols <- c(colnames(counts_all)[1:7], unique_vntg)
counts_all <- counts_all[keep_cols]
fpkm_all <- fpkm_all[keep_cols]
tpm_all <- tpm_all[keep_cols]

#-------------- Write data --------------#
write.csv(counts_all, "data/TMA36_project/RNA_Seq/processed/rawcounts_hg38.csv", row.names = FALSE)
write.csv(fpkm_all, "data/TMA36_project/RNA_Seq/processed/fpkm_hg38.csv", row.names = FALSE)
write.csv(tpm_all, "data/TMA36_project/RNA_Seq/processed/tpm_hg38.csv", row.names = FALSE)
