# Visualization

dendrogram_barplot <- function(data, dist_method = 'euclidean', 
    hclust_method = 'ward.D', corr_mat = FALSE, coul, xlabl="Cell type (% of sample)",
    ylabl="Patient ID", legend_p ='none', BPOrderAsDendrogram=T, bp_order,
    cex_names=0.75, cex_axis=0.75, cex_lab=1, cex_sub=1){
    
    # Dendrogram
    if(corr_mat){
        corr_pt <- Hmisc::rcorr(t(as.matrix(data)), type = 'spearman')
        dist_mat <- dist(corr_pt$r, method = dist_method)
    } else {
        dist_mat <- dist(data, method = dist_method)
    }
    hclust_avg <- hclust(dist_mat, method = hclust_method)
    par(mar=c(2,7,4,2), lwd=2, mgp=c(3,1,0), las=1, font.lab=2)
    plot(hclust_avg,cex = 0.8, hang = -1, main = paste0(hclust_method, ' linkage'))
    
    # Barplot
    #coul = coul
    if(BPOrderAsDendrogram){
      data <- data[hclust_avg$order,] #dendrogram order
    }else{
      data <- data[bp_order,]
    }
    
    data <- t(as.matrix(data))
    
    par(mgp=c(1.5,0.25,0), las=1, mar = c(3, 5, 1, 7.5), xpd=TRUE, font.lab=2, lwd=1) # axis label locations
    barplot(data,
            col=coul ,
            border='white',
            horiz=TRUE,
            cex.names = cex_names,
            cex.axis = cex_axis,
            cex.lab = cex_lab,
            #cex.sub = cex_sub,
            legend = FALSE)
    title(xlab=xlabl, mgp = c(1.5, 0.5, 0))
    title(ylab=ylabl,  mgp = c(3, 0.1, 0))
    # legend('topright', legend = rownames(data), fill = coul, ncol = 1, inset=c(-0.3,0),
    #            cex = 0.75) 

}

###############################################################################
# UMAP
###############################################################################

ClusterUMAP_plot <- function(data, 
                             density = F,
                             color_by_cluster = T, 
                             color_by_protein = F,
                             color_by_continuous = F, 
                             cluster_col, ft_cols, 
                             pallete = F,
                             plot_clusnames = F,
                             plot_colors,
                             plot_title,
                             lg_names,
                             x_lim = c(-10,10),
                             y_lim = c(-10,10),
                             nrow_prot = 2,
                             heights_prot = unit(c(1.3,1.3), c("in", "in")),
                             title_size_prot = 10, epi_clus=F){
  
  if(density){
    density_by=T
    test <- umap_smp[c('UMAP1', 'UMAP2', 'n_op2')]
    test$n_op2 <- as.character(test$n_op2)
    test$n_op2[test$n_op2=='ind'] = 'Indolent'
    test$n_op2[test$n_op2=='int'] = 'Intermediate'
    test$n_op2[test$n_op2=='agg'] = 'Aggressive'
    test$n_op2 <- factor(test$n_op2, levels = c('Indolent', 'Intermediate', 'Aggressive'))
    p_cl <- ggplot(test, aes(x = UMAP1, y = UMAP2)) + 
      geom_point(shape = 20, size = 0.01) + 
      ggtitle(plot_title) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line =  element_blank(),
            legend.position = "right",
            legend.title=element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            aspect.ratio = 1, axis.text = element_text(colour = 1, size = 10),
            legend.background = element_blank(),
            legend.spacing.y = unit(0, "mm"),
            legend.spacing.x = unit(0, "mm"),
            legend.box.background = element_rect(colour = "black"),
            legend.text = element_text(size = 8),
            plot.title = element_text(size=16, face="bold"),
            legend.key.size = unit(0.5, "cm"),
            legend.key.width = unit(0.5,"cm"),
            axis.title=element_text(size=10)) +
      #stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", n=100) +
      stat_density_2d(aes(fill = after_stat(piece)), geom = "polygon", n=100, bins=10, contour = T) +
      #xlim(x_lim) +
      #ylim(y_lim) +
      scale_x_continuous(limits = x_lim, expand = c(0, 0)) +
      scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
      scale_fill_viridis_c(option = 'magma') #https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html
    p_cl2 <- p_cl+facet_wrap(vars(n_op2),scales='free')
    p_cl_ls <- list(p_cl, p_cl2) #plot(p_cl), plot(p_cl2)
    return(p_cl_ls)
  }
  
  if(color_by_cluster){
    test <- data[c(cluster_col, 'UMAP1', 'UMAP2')]
    test[,cluster_col] <- as.factor(test[,cluster_col])
    
    p_cl <- ggplot(test) + 
      geom_point(aes_string(x='UMAP1', y='UMAP2', colour = test[,cluster_col]), 
                 shape = 20, size = 0.8, alpha = 0.4) +
      ggtitle(plot_title) +
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line =  element_blank(),
            legend.position = "right",
            legend.title=element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            aspect.ratio = 1, axis.text = element_text(colour = 1, size = 10),
            legend.background = element_blank(),
            legend.spacing.y = unit(0, "mm"),
            legend.spacing.x = unit(0, "mm"),
            legend.box.background = element_rect(colour = "black"),
            legend.text = element_text(size = 8),
            plot.title = element_text(size=16, face="bold"),
            legend.key.size = unit(0.5, "cm"),
            legend.key.width = unit(0.5,"cm"),
            axis.title=element_text(size=10)) +
      scale_x_continuous(limits = x_lim, expand = c(0, 0)) +
      scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
      guides(colour = guide_legend(override.aes = list(size=3, alpha=1, pch=15))) #override legend features
    
    if(pallete){
      #p_cl <- p_cl + scale_color_brewer(palette="Set3", labels = lg_names) + scale_color_hue(l=60, c=65)
      p_cl <- p_cl + scale_color_jcolors(palette = 'pal8', labels = lg_names)
    } else {
      p_cl <- p_cl + scale_color_manual(values=plot_colors, labels = lg_names)
    }
    
    if (plot_clusnames){
      edata <- cbind(test$UMAP1, test$UMAP2, cluster_col=test[cluster_col])
      colnames(edata) <- c('x', "y", "z")
      center <- aggregate(cbind(x,y) ~ z, data = edata, median)
      if(epi_clus){
        lb <- sapply(strsplit(as.character(center[,1]), "_"), "[[", 2)
      }else{
          lb <- center[,1]
        }
      p_cl <- p_cl + annotate("text", label = lb, 
                              x=center[,2], y = center[,3], size = 6, colour = "black", fontface = 'bold')
    }
    return(p_cl)
  }
  
  if(color_by_protein){
    data2 <- data[,ft_cols]
    data <- cbind(data2, data[c('UMAP1', 'UMAP2')])
    #data <- data[c(colnames(data)[ft_cols], 'UMAP1', 'UMAP2')]
    test <- melt(data, id= c('UMAP1', 'UMAP2'))
    var_list <- unique(test$variable)
    
    #x <- sapply(strsplit(var_list, "_"), "[[", 2)
    pl <- list()
    for(i in seq_along(var_list)) {
      p <- ggplot(subset(test,test$variable==var_list[i]), aes(x=UMAP1, y=UMAP2, colour = value)) + 
        geom_point(shape = 20, alpha=0.4, size = 0.3) + 
        theme_bw() + 
        ggtitle(sapply(strsplit(as.character(var_list[i]), "_"), "[[", 2)) + 
        scale_colour_gradient(low = "grey72", high = "blue") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.background = element_blank(), axis.line =  element_blank(),
              #axis.ticks=element_blank(), 
              axis.title=element_blank(),
              axis.text=element_blank(), aspect.ratio = 1,
              legend.position = "none",
              plot.title = element_text(size=title_size_prot)) +
      scale_x_continuous(limits = x_lim, expand = c(0, 0)) +
      scale_y_continuous(limits = y_lim, expand = c(0, 0))
      
      pl[[i]] <- p
    }
    umap_pl <- grid.arrange(grobs=pl, nrow=nrow_prot, heights=heights_prot)
    return(umap_pl)
  }

  if(color_by_continuous){
    test <- data[c(cluster_col, 'UMAP1', 'UMAP2')]
    
    p_cl <- ggplot(test) + 
      geom_point(aes_string(x='UMAP1', y='UMAP2', colour = test[,cluster_col]), 
                 shape = 20, size = 0.8, alpha = 0.4) +
      ggtitle(plot_title) +
      theme_bw()+
      scale_colour_gradient2(low = "#3498DB", mid= 'white', high = "#c75264", midpoint = 0.5) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line =  element_blank(),
            legend.position = "right",
            legend.title=element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            aspect.ratio = 1, axis.text = element_text(colour = 1, size = 10),
            legend.background = element_blank(),
            legend.spacing.y = unit(0, "mm"),
            legend.spacing.x = unit(0, "mm"),
            legend.box.background = element_rect(colour = "black"),
            legend.text = element_text(size = 8),
            plot.title = element_text(size=16, face="bold"),
            legend.key.size = unit(0.5, "cm"),
            legend.key.width = unit(0.5,"cm"),
            axis.title=element_text(size=10)) +
      scale_x_continuous(limits = x_lim, expand = c(0, 0)) +
      scale_y_continuous(limits = y_lim, expand = c(0, 0)) 
  
  return(p_cl)
  }  
}


###############################################################################
# BARPLOTS of subtypes
###############################################################################

ClassDistribution <- function(data, cond_col, celltype_col, cluster_col, celltype_nm, ptID_col){
    data <- data[which(data[,celltype_col] == celltype_nm), ]
    
    # Computing table of cluster %by behavior condition
    clust <- unique(data[,cluster_col])
    cond <- unique(data[,cond_col])
    prcnt_cond <- data.frame(matrix(nrow = length(clust), ncol = length(cond) , dimnames=list(clust, cond)))
    
    for(i in 1:nrow(prcnt_cond)){
      k <- which(data[,cluster_col] == clust[i])
      k <- table(data[k,cond_col])/length(k)
      k <- k[cond]
      k[is.na(k)] <- 0
      prcnt_cond[i,] <- k
    }

    # Computing table of cluster cell number by behavior condition
    num_clus <- data.frame(matrix(nrow = length(clust), ncol = 1 , dimnames=list(clust, 'cell_num')))
    for(i in 1:nrow(num_clus)){
      k <- length(which(data[,cluster_col] == clust[i]))
      num_clus[i,1] <- k
    }

    # Computing table of cluster by patient
    ptids <- unique(data[,ptID_col])
    prcnt_pt <- data.frame(matrix(nrow = length(clust), ncol = length(ptids) , dimnames=list(clust, ptids)))
    for(i in 1:nrow(prcnt_pt)){
      k <- which(data[,cluster_col] == clust[i])
      k <- table(data[k,ptID_col])/length(k)
      k <- k[ptids]
      k[is.na(k)] <- 0
      prcnt_pt[i,] <- k
    }

    res <- list("prcnt_cond"=prcnt_cond, "num_clus"=num_clus, "prcnt_pt"=prcnt_pt)

    # create a list with objetcs                  
    return(res)
}



bp_stclusters <- function(dt, celltype_nm_plot, plot_type){
    library(reshape2)
    library(ggplot2)
    library(RColorBrewer)
    if(!(plot_type %in% c('cond_frac', 'pt_frac', 'cell_num'))){
        stop('error')
    }

    dt['cluster'] <- sapply(strsplit(rownames(dt), "_"), "[[", 2)
    dt <- melt(dt, id.vars = c( "cluster"))
    dt$cluster <- factor(dt$cluster, levels = sort(as.numeric(unique(dt$cluster)), decreasing=T))


    if(plot_type == 'cond_frac'){
        if(length(unique(dt$variable))>2){
            dt$variable <- factor(dt$variable, levels = c('ind', 'int', 'agg'))
            col <- c("#3498DB", "grey72", "#c75264")
            lb <- c("Indolent", "Intermediate", "Aggressive")
        } else {
            dt$variable <- factor(dt$variable, levels = c('ind', 'agg'))
            col <- c("#3498DB", "#c75264")
            lb <- c("Indolent", "Aggressive")
        }
        p <- ggplot(dt, aes(fill=variable, y=value, x=cluster)) + 
                geom_bar(position="fill", stat="identity") +
                coord_flip() +
                scale_fill_manual(values=col, name = "Behavior", labels = lb) +
                xlab(celltype_nm_plot) + 
                ylab("Fraction of cells") +
                theme_minimal()+
                theme(legend.position="top")+
                guides(fill=guide_legend(title.position = "top", title.hjust=0.5))           
    }

    if(plot_type == 'pt_frac'){
        dt$variable <- as.character(sapply(strsplit( as.character(unique(dt$variable)), "X"), "[[", 2))
        set.seed(101)
        n <- length(unique(dt$variable))
        qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
        col <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
        lb <- unique(dt$variable)
        p <- ggplot(dt, aes(fill=variable, y=value, x=cluster)) + 
                geom_bar(position="fill", stat="identity") +
                coord_flip() +
                scale_fill_manual(values=col, name = "Patient", labels = lb) +
                xlab(celltype_nm_plot) + 
                ylab("Fraction of cells") +
                theme_minimal()+
                theme(legend.position="top")+
                guides(fill=guide_legend(ncol=10, title.position = "top", title.hjust=0.5))
    }

    if(plot_type == 'cell_num'){
        p <- ggplot(dt, aes(x=cluster, y=value)) +
                geom_bar(stat="identity", fill="#07315e")+
                theme_minimal() +
                coord_flip() +
                xlab(celltype_nm_plot) + 
                ylab("Number of cells")
    }
    return(p)

}

main_bp <- function(data, cond_col, celltype_col, cluster_col, celltype_nm, ptID_col, celltype_nm_plot){

    dt_ls <- ClassDistribution(
        data = data, 
        cond_col = cond_col, 
        celltype_col = celltype_col, 
        cluster_col = cluster_col, 
        celltype_nm = celltype_nm, 
        ptID_col = ptID_col)
    
    cond <- bp_stclusters(dt=dt_ls$prcnt_cond, celltype_nm_plot=celltype_nm_plot, plot_type='cond_frac')
    frac <- bp_stclusters(dt=dt_ls$prcnt_pt, celltype_nm_plot=celltype_nm_plot, plot_type='pt_frac')
    cnum <- bp_stclusters(dt=dt_ls$num_clus, celltype_nm_plot=celltype_nm_plot, plot_type='cell_num')

    bp_ls <- list(cond=cond, frac=frac, cnum=cnum)
    
    return(bp_ls)
}

frac_hm <- function(prcnt_dt, CDE, class_col, scale_v = TRUE){
    prcnt_dt[is.na(prcnt_dt)] <- 0
    lb <- rep('Indolent', nrow(prcnt_dt))
    lb[which(CDE[,class_col] == 'agg')] <- 'Aggressive'
    levels = c('Indolent', 'Aggressive')
    colr = list(Behavior = c('Indolent' = '#3498DB', 'Aggressive'= '#c75264'))

    if(length(unique(CDE[,class_col]))>2){
        lb[which(CDE[,class_col] == 'int')] <- 'Intermediate'
        levels = c('Indolent', 'Intermediate', 'Aggressive')
        colr = list(Behavior = c('Indolent' = '#3498DB', 'Intermediate' = 'grey72', 'Aggressive'= '#c75264'))
    }
    
    ha = HeatmapAnnotation(
        Behavior = factor(lb, levels = levels),
        col = colr,
        simple_anno_size = unit(0.5, "cm")
    )
    if(scale_v){
      dt <- t(as.matrix(scale(prcnt_dt)))
    }else{
      dt <- t(as.matrix(prcnt_dt))
      }
    
    dt[is.na(dt)] <- 0
    Heatmap(dt, name = "z-score",
            heatmap_legend_param = list(color_bar = "continuous"), 
            row_names_gp = gpar(fontsize = 8),
            column_names_gp = gpar(fontsize = 8),
            column_split = lb,
            top_annotation = ha)
}


frac_boxplot <- function(prcnt_dt, CDE, class_col, indvsagg=TRUE){
    prcnt_dt[is.na(prcnt_dt)] <- 0
    
    if(indvsagg){
      lb <- rep('Ind', nrow(prcnt_dt))
      lb[which(CDE[,class_col] == 'agg')] <- 'Agg'
      levels = c('Ind', 'Agg')
      comp =  list(c("Ind", "Agg"))
      colr = c("#3498DB", "#c75264")
      
      if(length(unique(CDE[,class_col]))>2){
        lb[which(CDE[,class_col] == 'int')] <- 'Int'
        levels = c('Ind', 'Int', 'Agg')
        comp =  list(c("Ind", "Agg"), c("Ind", "Int"), c("Int", "Agg"))
        colr = c("#3498DB", "grey72","#c75264")
        
      }
      
    }else{
      #CDE[,class_col] <- as.factor(CDE[,class_col])
      levels <- sort(unique(CDE[,class_col]))
      lb <- rep(levels[1], nrow(prcnt_dt))
      lb[which(CDE[,class_col] == levels[2])] <- levels[2]
      comp = list(c(levels))
      colr = c("#3498DB", "#c75264")
    }
    
    prcnt_dt['Group'] <- factor(lb, levels = levels)
    prcnt_dt['pt_ID'] <- CDE$pt_ID
    prcnt_dt <- reshape2::melt(prcnt_dt, id.vars = c("pt_ID", "Group"))
    colnames(prcnt_dt)[4] <- 'Fraction'
    
    sigFunc = function(x){
      if(x < 0.001){"***"} 
      else if(x < 0.01){"**"}
      else if(x < 0.05){"*"}
      else{NA}}
    
    ggplot(prcnt_dt, aes(x=Group, y=Fraction, color = Group)) +
        geom_boxplot() +
        #ylim(0,1)+
        ggsignif::geom_signif(comparisons = comp, 
                              map_signif_level= sigFunc,
                              margin_top = 0.05,color ='black',
                              step_increase = 0.06) +
        facet_wrap(~variable, scales='free_y') +
        scale_y_continuous(expand = expansion(mult = c(0.05, .2)))+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5, size=22),
              legend.position='none',
              #axis.text.x = element_text(angle = 45, hjust = 1)
              )+
        scale_color_manual(values=colr, name = "Group", labels = lb)

}


hm_median <- function(data){
    source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/20_ClustAnnot_functions.R")
    data$cluster <- as.numeric(gsub(".*_",  "",data$cluster))
    data <- ClusterEval_data(data[,1:(ncol(data)-1)], eval_type = 'heatmap')
    colnames(data) <- gsub(".*_",  "",colnames(data))

    library(gplots)
    scale_max = max(data)
    heat_palette_med <- colorRampPalette(c("black", "yellow","#FAF7C9"))
    pairs.breaks_med <- c(seq(0, scale_max/6.6, by = 0.1),
                          seq(scale_max/6.6, scale_max/3.3, by = 0.1),
                          seq(scale_max/3.3, scale_max, by = 0.1))

    heatmap.2(data,
              main = "Median protein expression",
              dendrogram = "both",
              Rowv = TRUE,
              Colv = TRUE,
              breaks = pairs.breaks_med,
              revC = FALSE,
              symkey = FALSE,
              symbreaks = FALSE,
              scale = "none",
              cexRow = 1.2,
              cexCol = 1.2,
              key = TRUE,
              col = heat_palette_med,
              trace = "none",
              density.info = 'none',
              sepcolor="#424242",
              margins = c(7,3),
              colsep=1:ncol(data),
              rowsep=1:nrow(data),
              sepwidth=c(0.005,0.005),
              keysize = 1.2,
              key.title = 'Intensity',
              key.xlab= "Arcsinh Transform",
              extrafun = box(lty = "solid"),
              srtCol=90,
              lhei=c(1,4), lwid=c(5,25))

}

protein_corr <- function(corr_ls) {
    #par(cex=0.9) # size
    #par(lwd = 1.2) # line width
    r <- corr_ls$r
    pv <- corr_ls$P
    colnames(r) <- gsub(".*_",  "",colnames(r))
    rownames(r) <- colnames(r)
    colnames(pv) <- gsub(".*_",  "",colnames(pv))
    rownames(pv) <- colnames(pv)
    pv[is.na(pv)] = 0
    corrplot::corrplot(r, type="upper", order='hclust', tl.col = "black", 
        tl.srt = 45, p.mat = pv , sig.level = 0.05, 
        insig = "blank", method = 'color', addCoef.col="black", 
        number.font = 1, number.cex = 0.8, addgrid.col = 'grey', tl.cex = 1)  
}

hist_prot <- function(data){
    data <- data[,1:(ncol(data)-1)]
    colnames(data) <- gsub(".*_",  "",colnames(data))
    data <- reshape2::melt(data,  id.vars = c('cluster'))
    data$cluster <- as.factor(data$cluster)
    ggplot(data, aes(x=value, group=cluster, color=cluster)) +
        geom_density(adjust=1.5) + facet_wrap( ~ variable, scales="free")
}
