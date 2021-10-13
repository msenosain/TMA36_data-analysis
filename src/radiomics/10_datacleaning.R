###############################################################################
# Radiomics
###############################################################################
## 1. Healthmyne
m_HM <- read.csv('data/TMA36_project/Radiomics/processed/TMA36_HealthMyne_Khushbu.csv')
rownames(m_HM) <- m_HM[,1]
m_HM <- m_HM[-which(is.na(m_HM[,2])),]
m_HM <- m_HM[,15:ncol(m_HM)]
nona <- colnames(m_HM)[apply(m_HM, 2, anyNA)]
m_HM <- m_HM[,-which(colnames(m_HM) %in% nona)]

# Center and scale
#m_HM <- scale(m_HM)

write.csv(m_HM, 'data/TMA36_project/Radiomics/processed/rad_healthmyne.csv')
