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

CDE['Death_st'] <- 'No'
k1 <- which(is.na(CDE$Death_Date)==FALSE)
CDE[k1,'Death_st'] <- 'Yes'

CDE['Recurrence_st'] <- 'No'
k2 <- which(is.na(CDE$Recurrence_Date)==FALSE)
CDE[k2,'Recurrence_st'] <- 'Yes'

CDE['Progression_st'] <- 'No'
k3 <- which(is.na(CDE$Progression_Date)==FALSE)
CDE[k3,'Progression_st'] <- 'Yes'

CDE['DRP_st'] <- 'No'
CDE[unique(c(k1,k2,k3)),'DRP_st'] <- 'Yes'

write.csv(CDE, 'data/TMA36_project/CDE/CDE_TMA36_2021SEPT21_DR_MF.csv', row.names = FALSE)
