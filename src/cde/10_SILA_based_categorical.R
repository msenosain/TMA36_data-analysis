CDE <- read.csv('data/TMA36_project/CDE/CDE_TMA36_2021SEPT21_DR.csv')

ind <- which(CDE$SILA <= 0.4)
int <- which(CDE$SILA > 0.4 & CDE$SILA <= 0.6)
agg <- which(CDE$SILA > 0.6)
CDE['n_op2'] <- c(rep(0, nrow(CDE)))
CDE[ind, 'n_op2'] <- 'ind'
CDE[int, 'n_op2'] <- 'int'
CDE[agg, 'n_op2'] <- 'agg'

colnames(CDE)[1]<- 'pt_ID'

date_cols <- grep('Date|date|dob|LDKA', colnames(CDE))
CDE[,date_cols] <- apply(CDE[,date_cols], 2, function(x) format(strptime(as.character(x), "%m/%d/%Y"), "%Y-%m-%d"))

write.csv(CDE, 'data/TMA36_project/CDE/CDE_TMA36_2021SEPT21_DR_MF.csv', row.names = FALSE)
