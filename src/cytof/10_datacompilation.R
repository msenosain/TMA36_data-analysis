
######----------DEFINE FUNCTIONS----------######

datacomp <- function(CD90=TRUE, filename = filename, files_path){
    # Set working directory
    original_wd <- getwd()
    setwd(files_path)
    
    # Read files list
    fcs <- '.fcs|.FCS'

    files_list <- list.files(pattern=fcs)
    print('File list successfully created!')

    # Remove files of controls
    files_list <- files_list[!grepl("._0_.", files_list)]

    # Remove files of exluded patients
    files_list <- files_list[!grepl(".X.", files_list)]


    # Extract CyTOF dates
    cytof_dates <- as.Date(gsub("_.*$",  "",files_list), "%y%m%d")

    # Extract patient ID
    #pt_ID <- as.numeric(sapply(strsplit(files_list, "_"), "[[", 3))
    pt_ID <- sapply(strsplit(files_list, "_"), "[[", 3)

    CANARY <- sapply(strsplit(files_list, "_"), "[[", 4)

    # State if it has CD90 antibody or not
    if(CD90){
        CD90 <- rep('Yes', length(files_list)) # add a promt to choose with or without
    } else {
        CD90 <- rep('No', length(files_list)) # add a promt to choose with or without
    }
    

    # Creat a dataframe with all vectors
    ref <- cbind.data.frame('pt_ID' = pt_ID, 'CANARY' = CANARY, 
        'CyTOF_date' = as.Date(cytof_dates), "CD90" = CD90)
    #ref$CyTOF_date <- as.Date(ref$CyTOF_date)


    # Read first file to correct column order later
    smp <- flowCore::read.FCS(files_list[1], transformation = FALSE)
    descrp <- smp@parameters@data$desc
    smp <- data.frame(smp@exprs)
    col_nms <- colnames(smp)

    big_df <- data.frame(matrix(ncol = ncol(smp)+3, nrow=0))
    colnames(big_df) <- c(col_nms, 'pt_ID', 'CANARY', 'CD90')

    for (i in 1:length(files_list)){
        # Read exprs from FCS file
        dt <- data.frame(flowCore::read.FCS(files_list[i], transformation = FALSE)@exprs)

        # Order columns
        dt <- dt[col_nms]

        dt['pt_ID'] <- rep(pt_ID[i], nrow(dt))
        dt['CANARY'] <- rep(CANARY[i], nrow(dt))
        dt['CD90'] <- rep(CD90[i], nrow(dt))

        big_df <- rbind(big_df, dt)
    }

    colnames(big_df)<- c(descrp, 'pt_ID', 'CANARY', 'CD90')

    # Save RData file
    save(big_df, ref, file = paste0(filename, '.RData'))
    
    # Reset working directory
    setwd(dir = original_wd)
}


datacomp_ctl <- function(filename = filename, files_path){
    # Set working directory
    original_wd <- getwd()
    setwd(files_path)
    
    # Read files list
    fcs <- '.fcs|.FCS'

    files_list <- list.files(pattern=fcs)
    print('File list successfully created!')

    # Remove files of controls
    files_list <- files_list[grepl("._0_.", files_list)]


    # Extract CyTOF dates
    cytof_dates <- as.Date(gsub("_.*$",  "",files_list), "%y%m%d")
    

    # Creat a dataframe with all vectors
    ref <- cbind.data.frame('batch_ID' = c(1:length(files_list)), 
        'CyTOF_date' = as.Date(cytof_dates))


    # Read first file to correct column order later
    smp <- flowCore::read.FCS(files_list[1], transformation = FALSE)
    descrp <- smp@parameters@data$desc
    smp <- data.frame(smp@exprs)
    col_nms <- colnames(smp)

    big_df <- data.frame(matrix(ncol = ncol(smp)+1, nrow=0))
    colnames(big_df) <- c(col_nms, 'CyTOF_date')

    for (i in 1:length(files_list)){
        # Read exprs from FCS file
        dt <- data.frame(flowCore::read.FCS(files_list[i], transformation = FALSE)@exprs)

        # Order columns
        dt <- dt[col_nms]

        dt['CyTOF_date'] <- rep(cytof_dates[i], nrow(dt))

        big_df <- rbind(big_df, dt)
    }

    colnames(big_df)<- c(descrp, 'CyTOF_date')

    # Save RData file
    save(big_df, ref, file = paste0(filename, '.RData'))
    
    # Reset working directory
    setwd(dir = original_wd)
}

######----------RUN FUNCTIONS----------######

# Samples with CD90
datacomp(CD90=TRUE, filename = 'withCD90', 
         files_path = 'data/TMA36_project/CyTOF/processed/all/renamed/normed/normed_samples/withCD90/output/')

# Samples without CD90
datacomp(CD90=FALSE, filename = 'woCD90', 
         files_path = 'data/TMA36_project/CyTOF/processed/all/renamed/normed/normed_samples/withoutCD90/output/')
    
# Controls
datacomp_ctl(filename = 'CyTOF_ADC_controls', files_path = 'data/TMA36_project/CyTOF/processed/all/renamed/normed/normed_samples/both_output/')

######----------MOVE FILES----------######
dir.create("data/TMA36_project/CyTOF/processed/Data_paper2/withCD90", recursive = TRUE)
dir.create("data/TMA36_project/CyTOF/processed/Data_paper2/woCD90")
dir.create("data/TMA36_project/CyTOF/processed/Data_paper2/controls")

filesstrings::file.move("data/TMA36_project/CyTOF/processed/all/renamed/normed/normed_samples/withCD90/output/withCD90.RData",
          "data/TMA36_project/CyTOF/processed/Data_paper2/withCD90")
filesstrings::file.move("data/TMA36_project/CyTOF/processed/all/renamed/normed/normed_samples/withoutCD90/output/woCD90.RData",
                        "data/TMA36_project/CyTOF/processed/Data_paper2/woCD90")
filesstrings::file.move("data/TMA36_project/CyTOF/processed/all/renamed/normed/normed_samples/both_output/CyTOF_ADC_controls.RData",
                        "data/TMA36_project/CyTOF/processed/Data_paper2/controls")


