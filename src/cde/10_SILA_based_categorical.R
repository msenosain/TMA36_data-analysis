CDE <- read.csv('data/TMA36_project/CDE/CDE_TMA36_2020FEB25_SA.csv')

ind <- which(CDE$SILA <= 0.4)
int <- which(CDE$SILA > 0.4 & CDE$SILA <= 0.6)
agg <- which(CDE$SILA > 0.6)
CDE['n_op2'] <- c(rep(0, nrow(CDE)))
CDE[ind, 'n_op2'] <- 'ind'
CDE[int, 'n_op2'] <- 'int'
CDE[agg, 'n_op2'] <- 'agg'

write.csv(CDE, 'data/TMA36_project/CDE/CDE_TMA36_2020FEB25_SA_MF.csv', row.names = FALSE)
