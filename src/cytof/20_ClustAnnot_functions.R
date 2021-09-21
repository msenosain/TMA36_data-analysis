# Select columns of interest
# Transform data
# Clustering
# plot clusters (HM, tSNE)
# Cluster annotation
# add new column with cell types to big df

# Subtype clustering

clustering <- function(data, n_clusters = 10, iterations = 200, seed = 101){
    
    # Subseting data for clustering and transforming data
    print(as.matrix(colnames(data)))
    prompt <- "Enter the column INDICES of the training features (separated by single space only, no comas allowed) \n"
    features <- as.numeric(strsplit(readline(prompt), " ")[[1]])
    features <- colnames(data)[features]
    df_cluster <- denoisingCTF::t_asinh(data[features])

    # Clustering
    set.seed(seed)
    cl_data <- stats::kmeans(df_cluster, centers = n_clusters, iter.max = iterations)$cluster
    df_cluster <- cbind(df_cluster, 'cluster'= cl_data)

    return(df_cluster)
}




ClusterEval_data <- function(df_cluster, eval_type = c('heatmap', 'tSNE', 'UMAP'), 
    sample_size = 50000){
    # Heatmap
    if(eval_type == 'heatmap'){
        cl_median <- aggregate(df_cluster, list(df_cluster$cluster), median)[,-1]
        cl_median <- round(as.matrix(cl_median[,1:ncol(cl_median)-1]), digits = 3)
        return(cl_median)
    }
    
    # tSNE
    if(eval_type == 'tSNE'){
        smpl <- dplyr::sample_n(df_cluster, size=sample_size, replace=F)
        tsne_smp <- Rtsne::Rtsne(smpl[,-ncol(df_cluster)],  check_duplicates = FALSE)

        # Editing df
        tsne_smp <- as.data.frame(cbind(smpl,
                                         tSNE1 = tsne_smp$Y[,1],
                                         tSNE2 = tsne_smp$Y[,2]))
        return(tsne_smp)
    }

    # U-MAP
    if(eval_type == 'UMAP'){
        smpl <- dplyr::sample_n(df_cluster, size=sample_size, replace=F)
        umap_smp <- umap::umap(smpl[,-ncol(df_cluster)])

        # Editing df
        umap_smp <- as.data.frame(cbind(smpl,
                                         UMAP1 = umap_smp$layout[,1],
                                         UMAP2 = umap_smp$layout[,2]))
        return(umap_smp)
    }
}

ClusterEval_plot <- function(data, data_type = c('heatmap', 'tSNE', 'UMAP')){
    # Heatmap
    if(data_type == 'heatmap'){
        library(gplots)
        scale_max = max(data)
        heat_palette_med <- colorRampPalette(c("black", "yellow","#FAF7C9"))
        pairs.breaks_med <- c(seq(0, scale_max/6.6, by = 0.1),
                              seq(scale_max/6.6, scale_max/3.3, by = 0.1),
                              seq(scale_max/3.3, scale_max, by = 0.1))

        hm_pl <- heatmap.2(data,
                  main = "Median protein expression",
                  dendrogram = "both",
                  Rowv = TRUE,
                  Colv = TRUE,
                  breaks = pairs.breaks_med,
                  revC = FALSE,
                  symkey = FALSE,
                  symbreaks = FALSE,
                  scale = "none",
                  cexRow = 0.8,
                  cexCol = 0.8,
                  key = TRUE,
                  col = heat_palette_med,
                  trace = "none",
                  density.info = 'none',
                  sepcolor="#424242",
                  margins = c(6,10),
                  colsep=1:ncol(data),
                  rowsep=1:nrow(data),
                  sepwidth=c(0.005,0.005),
                  keysize = 1,
                  key.title = 'Intensity',
                  key.xlab= "Arcsinh Transform",
                  extrafun = box(lty = "solid"),
                  cellnote = data,
                  notecol = 'red',
                  srtCol=45
                  )
        return(hm_pl)
    }

    if(data_type == 'tSNE'){
        # Plotting
        library(ggplot2)
        library(gridExtra)
        library(reshape2)
        test <- melt(data, id= c('cluster', 'tSNE1', 'tSNE2'))
        var_list <- unique(test$variable)
        pl <- list()

        # Defining cluster centers
        edata <- cbind(data$tSNE1, data$tSNE2, data$cluster)
        colnames(edata) <- c('x', "y", "z")
        center <- aggregate(cbind(x,y) ~ z, data = edata, median)

        p_cl <- ggplot(data, aes(x=tSNE1, y=tSNE2, colour = factor(cluster)))+
          geom_point(alpha=0.3) + theme_bw() + ggtitle('cluster')+
          annotate("text", label = center[,1], x=center[,2], y = center[,3],
                   size = 6, colour = "black", fontface = 'bold')
        pl[[1]] <- p_cl
        for(i in seq_along(var_list)) {
          p <- ggplot(subset(test,test$variable==var_list[i]), aes(x=tSNE1, y=tSNE2, colour = value)) +
            geom_point(alpha=0.3) + theme_bw() + ggtitle(var_list[i]) +
            scale_colour_gradient(low = "gray75", high = "blue")
          pl[[i+1]] <- p
        }

        tsne_pl <- grid.arrange(grobs=pl)

        return(tsne_pl)
    }

    if(data_type == 'UMAP'){
        # Plotting
        library(ggplot2)
        library(gridExtra)
        library(reshape2)
        test <- melt(data, id= c('cluster', 'UMAP1', 'UMAP2'))
        var_list <- unique(test$variable)
        pl <- list()

        # Defining cluster centers
        edata <- cbind(data$UMAP1, data$UMAP2, data$cluster)
        colnames(edata) <- c('x', "y", "z")
        center <- aggregate(cbind(x,y) ~ z, data = edata, median)

        p_cl <- ggplot(data, aes(x=UMAP1, y=UMAP2, colour = factor(cluster)))+
          geom_point(alpha=0.3) + theme_bw() + ggtitle('cluster')+
          annotate("text", label = center[,1], x=center[,2], y = center[,3],
                   size = 6, colour = "black", fontface = 'bold')
        pl[[1]] <- p_cl
        for(i in seq_along(var_list)) {
          p <- ggplot(subset(test,test$variable==var_list[i]), aes(x=UMAP1, y=UMAP2, colour = value)) +
            geom_point(alpha=0.3) + theme_bw() + ggtitle(var_list[i]) +
            scale_colour_gradient(low = "gray75", high = "blue")
          pl[[i+1]] <- p
        }

        umap_pl <- grid.arrange(grobs=pl)

        return(umap_pl)
    }
}

ClusterAnnotation <- function(data, df_cluster, ls_annotation, 
    annotation_col = 'annotation_col', cl_delete = FALSE, cl_delete_name = 'delete'){
    
    # Cluster annotation and merge
    ct_k <- data.frame(cbind(cluster=df_cluster$cluster, annotation=rep(NA, nrow(df_cluster))))
    for (i in 1:length(ls_annotation)){
      k <- which(df_cluster$cluster %in% ls_annotation[[i]])
      ct_k[k,'annotation'] <- names(ls_annotation)[i]
    }
    colnms <- colnames(data)
    data <- cbind(data, ct_k$annotation)
    colnames(data) <- c(colnms, annotation_col)

    if(cl_delete){
        # Remove clusters
        k <- which(data[,annotation_col] == cl_delete_name)
        data <- data[-k,]
    }

    return(data)
}



ClassAbundanceByPt <- function(data, ptID_col = 'pt_ID', class_col){
    ptids <- unique(data[,ptID_col])
    ct <- unique(data[,class_col])
    prcnt <- data.frame(matrix(ncol = length(ct), nrow = length(ptids), dimnames=list(ptids, ct)))

    for(i in 1:nrow(prcnt)){
      k <- which(data[,ptID_col] == ptids[i])
      k <- table(data[k,class_col])/length(k)
      k <- k[ct]
      k[is.na(k)] <- 0
      prcnt[i,] <- k
    }

    return(prcnt)
}


findElbow <- function(data, max){
  n <- length(data)    
  data <- as.data.frame(cbind(1:n,data))
  colnames(data) <- c("X","Y")
  #r_plt <- c() #MF
  min_r <- Inf
  optimal <- 1
  for(i in 2:(n-1)){
    f1 <- stats::lm(Y~X,data[1:(i-1),])
    f2 <- stats::lm(Y~X,data[i:n,])
    r <- sum(abs(c(f1$residuals,f2$residuals)))
    #r_plt <- c(r_plt, r) #MF
    if(r < min_r){
      min_r <- r
      optimal <-i
    }
  }
  #plot(2:(max-1), r_plt, type="b", xlab="Number of Clusters", 
       #ylab="Residuals") #MF
  #abline(v = optimal, col = 'red')
  optimal
}


DetermineNumberOfClusters <- function(data,k_max,plot=FALSE,smooth=0.2,
                                      iter.max=50, seed = 101, ask_ft = T,
                                      arcsn_tr = T,
                                      ...){

  # Subseting data for clustering and transforming data
  if(ask_ft){
      print(as.matrix(colnames(data)))
      prompt <- "Enter the column INDICES of the training features (separated by single space only, no comas allowed) \n"
      features <- as.numeric(strsplit(readline(prompt), " ")[[1]])
      features <- colnames(data)[features]
      data <- data[features]
  } 

  if(arcsn_tr){
    data <- denoisingCTF::t_asinh(data)
  }
  # Clustering
  set.seed(seed)

  res <- sapply(1:k_max, 
              function(k){kmeans(data, k,iter.max = iter.max )$tot.withinss})
  
  for(i in 2:(k_max-1)){
    res[i] <- (1-smooth)*res[i]+(smooth/2)*res[i-1]+(smooth/2)*res[i+1]
  }
  
  
  if(plot) plot(1:k_max, res, type="b", xlab="Number of Clusters", 
                ylab="Within-cluster sum of squares", las=1)
  
  elbow <- findElbow(res,k_max)

  #abline(v = elbow, col = 'red')

  print(elbow)
}


# Plots

# Cluster Dendrogram

dendrogram_barplot <- function(data, dist_method = 'euclidean', 
    hclust_method = 'ward.D'){

    # Dendrogram
    dist_mat <- dist(data, method = dist_method)
    hclust_avg <- hclust(dist_mat, method = hclust_method)
    # par(mar=c(2,7,4,2), lwd=2)
    # plot(hclust_avg,cex = 0.8, hang = -1)

    if(ncol(data)>8){
        ncol_bp <- 6
        cex_bp <- 0.6
    } else {
        ncol_bp <- 4
        cex_bp <- 0.75
    }

    # Barplot
    coul = brewer.pal(9, "Set1")
    data <- data[hclust_avg$order,] #dendrogram order
    data <- t(as.matrix(data))
    par(las=1) # orientation, 1=horizontal
    # par(c(3.5,6,1,2)) #, lwd = 0.1) # mar: margins, lwd: line width
    # par(mgp=c(3,0.5,0)) # axis label locations
    barplot(data,
            col=coul ,
            border='white',
            horiz=TRUE,
            cex.names=1,
            cex.axis = 1,
            ylim = c(0,98))
    title(ylab="Patient ID", mgp=c(3.8,8,1), cex.lab=1.2)
    title(xlab="Cell type (% of sample)",  mgp=c(2,1,0), cex.lab=1.2)
    legend('top', legend = rownames(data), fill = coul, ncol = ncol_bp,
           cex = cex_bp)

}

corr_plot <- function(data, rcorr_type = 'spearman', p.adjust_method = 'BH'){

    res <- Hmisc::rcorr(as.matrix(data), type = rcorr_type) #for corr plot
    # corrplot
    corrected_pvals <- p.adjust(res$P, method = p.adjust_method)
    corrected_pvals <- matrix(corrected_pvals, nrow = ncol(res$P), 
        ncol = ncol(res$P))
    colnames(corrected_pvals)<- colnames(res$P)
    rownames(corrected_pvals)<- rownames(res$P)

    col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

    par(font = 2) # bold axis labels
    par(cex=0.9) # size
    par(lwd = 1.2) # line width
    corrplot::corrplot(res$r, type="upper", order='hclust', tl.col = "black", 
        tl.srt = 45, p.mat = corrected_pvals , sig.level = 0.05, 
        insig = "blank", method = 'color', addCoef.col="black", col=col(200),
        number.font = 1, number.cex = 0.8, addgrid.col = 'grey', tl.cex = 1)    
}



boxplots_m <- function(sbst, col_name, cell_pop){
    sbst <- sbst[which(sbst[,col_name] %in% cell_pop),]
    sbst <- subset(sbst, CANARY %in% c('G', 'P'))
    sbst <- denoisingCTF::t_asinh(sbst)

    markers <- colnames(sbst)[c(15, 17:31, 33:35, 37:48, 50, 51)]
    markers <- sapply(strsplit(markers, "_"), "[[", 2)
    markers <- gsub('-', '', markers)
    colnames(sbst)[c(15, 17:31, 33:35, 37:48, 50, 51)] <- markers

    library(dplyr)
    x <- sbst%>%
      group_by(CANARY)%>% 
      summarise_at(markers, median)
    x <- x[,-1]
    x <- t(x)
    k <- which(x[,1]<1.5 & x[,2]<1.5)

    markers <- markers[-k]

    sbst <- reshape2::melt(sbst)
    sbst <- subset(sbst, variable %in% markers)

    ggplot(sbst, aes(x=CANARY, y=value, color = CANARY)) +
       geom_boxplot() +
       ylim(0,12)+
       ggsignif::geom_signif(comparisons = list(c("G", "P")), 
              map_signif_level=TRUE) +
       facet_wrap(~variable) +
       ggtitle(cell_pop) +
       theme(plot.title = element_text(hjust = 0.5, size=22))

    # plotly::ggplotly(p) only if not using ggsignif pckg

}

edit_names <- function(df, col_name, var_oldname, var_newname){
  df[,col_name]<- as.character(df[,col_name])
  k <- which(df[,col_name] == var_oldname)
  df[k,col_name] <- var_newname
  df[,col_name]<- as.factor(df[,col_name])
  df
}

change_colname <- function(df, col_name, new_colname){
    i <- grep(col_name, colnames(df))
    cnms <- colnames(df)
    cnms[i] <- new_colname
    colnames(df) <- cnms
    df
}





