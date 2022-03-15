environment_set <- function(){
  library(EnhancedVolcano)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(DESeq2)
  library(sva)
  library(limma)
  library(ComplexHeatmap)
  library(fgsea)
  library(tibble)
  library(gage)
  library(ggplot2)
  library(org.Hs.eg.db)
  #library(clusterProfiler)
  library(survival)
  library(ggpubr)
  library(survminer)
  library(survRM2)
  # Check BiocManager::valid()
}

preprocess_rna <- function(raw_counts, batch_info, CDE,
    correct_batch=TRUE, 
    lowvargenesrm = TRUE,
    xychr_rm = TRUE,
    prot_coding_only = FALSE){

    # Remove low variance genes from counts
    if(lowvargenesrm){
        x<- raw_counts[,8:ncol(raw_counts)]
        idx <- edgeR::filterByExpr(x)
        raw_counts <- raw_counts[idx, ]
    }

    if(xychr_rm){
        idx <- which(raw_counts$Feature_chr %in% c('chrX', 'chrY'))
        raw_counts <- raw_counts[-idx, ]
    }

    if(prot_coding_only){
        idx <- which(raw_counts$Feature_gene_biotype == 'protein_coding')
        raw_counts <- raw_counts[idx, ]
    }    
    
    counts_only <- raw_counts[, 8:ncol(raw_counts)]
    rownames(counts_only) <- raw_counts$Feature

    # Normalize
    batch_info$Batch <- as.factor(batch_info$Batch)
    dds <- DESeqDataSetFromMatrix(
      countData = counts_only,
      colData = batch_info,
      design = ~Batch, tidy = F  
    )
    dds <- estimateSizeFactors(dds)
    sizeFactors(dds)
    vsd <- vst(dds)
    vsd_mat <- assay(vsd)

    # Assessing batch effect
    pbatch_bf <- plotPCA(vsd, "Batch") + labs(fill = "Batch") + ggtitle("Batch raw")
    
    # Matching with  CDE
    CDE <- CDE[match(batch_info$pt_ID, CDE$pt_ID),] # Select only pts with RNA Seq data

    # Save pre-processed data in a list
    dge_preprocessed <- list(batch_info=batch_info, raw_counts=raw_counts, 
                             counts_only=counts_only, vsd_mat=vsd_mat, CDE = CDE,
                             pbatch_bf=pbatch_bf)

    # Remove batch effect
    if(correct_batch){
        assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch=vsd$Batch)
        pbatch_af <- plotPCA(vsd, "Batch") + labs(fill = "Batch") + ggtitle("Batch after BE removal")
        vsd_mat <- assay(vsd)
        dge_preprocessed[['vsd_mat']] <- vsd_mat
        dge_preprocessed[['pbatch_af']] <- pbatch_af
    }
    
    dge_preprocessed
}

label_by_gene <- function(normalized_data, gene, metadata, extremes_only=FALSE){
    if(extremes_only){
        data_t <- data.frame(t(normalized_data))
        first_q <- quantile(data_t[,gene])[["25%"]]
        third_q <- quantile(data_t[,gene])[["75%"]]
        data_t$Condition <- "normal"
        data_t$Condition[data_t[,gene] < first_q] <- "low"
        data_t$Condition[data_t[,gene] > third_q] <- "high"

    } else {
        data_t <- data.frame(t(normalized_data))
        mid_q <- quantile(data_t[,gene])[["50%"]]
        data_t$Condition <- "low"
        data_t$Condition[data_t[,gene]  > mid_q] <- "high"

    }

    meta_data <- cbind(metadata, Condition=data_t$Condition)
    meta_data$Batch <- factor(meta_data$Batch)
    #meta_data$Gender <- factor(meta_data$Gender)
    meta_data$Condition <- factor(meta_data$Condition)

    meta_data
}


DE_analysis <- function(ls_preprocessed, 
    GeneBased=TRUE, 
    pDataBased=FALSE,
    NewCondition=FALSE,
    NewCondition_df, # condition label into batch_info
    cond_nm, # gene or column name
    two_levels=c('high','low'), # high low, dead alive...
    reference, # low, alive
    extremes_only=TRUE
    ){

    # Unlist
    batch_info = ls_preprocessed$batch_info
    raw_counts = ls_preprocessed$raw_counts
    CDE = ls_preprocessed$CDE
    counts_only = ls_preprocessed$counts_only
    vsd_mat = ls_preprocessed$vsd_mat

    message('Unlist done')


    if(GeneBased){
        meta_data <- label_by_gene(vsd_mat, gene=cond_nm, batch_info, extremes_only=extremes_only)
        meta_data <- meta_data %>% filter(Condition != "normal")
    }
    
    if(pDataBased){
        meta_data <- CDE %>% dplyr::select(pt_ID, all_of(cond_nm))
        meta_data$pt_ID <- as.character(meta_data$pt_ID)
        batch_info$pt_ID <- as.character(batch_info$pt_ID)
        meta_data <- inner_join(batch_info, meta_data, by='pt_ID')
        names(meta_data)[names(meta_data) == cond_nm] <- 'Condition'
        meta_data <- meta_data %>% filter(Condition %in% two_levels)
    }

    if(NewCondition){
        meta_data <- NewCondition_df
        print(meta_data)
        names(meta_data)[names(meta_data) == cond_nm] <- 'Condition'
        print(meta_data)
        meta_data$Condition <- as.character(meta_data$Condition)
        print(meta_data)
        meta_data <- meta_data %>% filter(Condition %in% two_levels)
        print(meta_data)
    }

    message('Labeling done')

    counts_only <- as.matrix(counts_only[,meta_data$Vantage_ID])
    CDE <- CDE %>% filter(pt_ID %in% meta_data$pt_ID)
    vsd_mat <- vsd_mat[,meta_data$Vantage_ID]

    message('Filtering done')

    dds <- DESeqDataSetFromMatrix(
      countData = counts_only,
      colData = meta_data,
      design = ~Condition, tidy = F) #design = ~Batch + Condition, tidy = F)

    message('Design done')

    dds <- estimateSizeFactors(dds)

    # Generate a transformed matrix with gene symbols
    ens2symbol <- data.frame(cbind(ENSEMBL=as.character(raw_counts$Feature), 
        symbol=as.character(raw_counts$Feature_gene_name)))
    vsd_mat_sym <- data.frame(vsd_mat) %>% mutate(gene=rownames(.)) %>% as_tibble()
    vsd_mat_sym <- inner_join(vsd_mat_sym, ens2symbol, by=c("gene"="ENSEMBL"))
    rownames(vsd_mat_sym) <- make.names(vsd_mat_sym$symbol, unique=TRUE)

    message('vsd symbols done')

    dds$Condition <- relevel(dds$Condition, ref = reference)
    dds <- DESeq(dds, parallel = F)

    message('DESeq done')

    res <- results(dds, contrast = c('Condition', two_levels))
    #res <- lfcShrink(dds, contrast = c('Condition',two_levels), res=res, type = 'normal')
    res_df <- data.frame(res) %>% mutate(gene=rownames(.)) %>% as_tibble()

    res_df <- res_df %>% 
      inner_join(., ens2symbol, by=c("gene"="ENSEMBL")) %>%
      data.frame()
    rownames(res_df) <- make.names(res_df$symbol, unique=TRUE)

    message('res symbols done')

    DE_res <- list(dds=dds, vsd_mat=vsd_mat, vsd_mat_sym=vsd_mat_sym, 
        ens2symbol=ens2symbol, res=res, res_df=res_df, 
        CDE=CDE, meta_data=meta_data)

    message('list done')

    DE_res
}

fgsea_analysis <- function(DE_res){

    # Unlist
    res_df=DE_res$res_df

    ranks <- res_df %>%
      dplyr::select(symbol,stat) %>%
      na.omit() %>% 
      distinct() %>% 
      group_by(symbol) %>% 
      summarize(stat=mean(stat)) %>%
      arrange(stat)

    ranks <- deframe(ranks)

    hallmark_path <- 'data/TMA36_project/RNA_Seq/GSEA/msigdb_v7.1_GMTs/h.all.v7.1.symbols.gmt'
    c1_path <- 'data/TMA36_project/RNA_Seq/GSEA/msigdb_v7.1_GMTs/c1.all.v7.1.symbols.gmt'
    c2_path <- 'data/TMA36_project/RNA_Seq/GSEA/msigdb_v7.1_GMTs/c2.all.v7.1.symbols.gmt'
    c3_path <- 'data/TMA36_project/RNA_Seq/GSEA/msigdb_v7.1_GMTs/c3.all.v7.1.symbols.gmt'
    c4_path <- 'data/TMA36_project/RNA_Seq/GSEA/msigdb_v7.1_GMTs/c4.all.v7.1.symbols.gmt'
    c5_path <- 'data/TMA36_project/RNA_Seq/GSEA/msigdb_v7.1_GMTs/c5.all.v7.1.symbols.gmt'
    c6_path <- 'data/TMA36_project/RNA_Seq/GSEA/msigdb_v7.1_GMTs/c6.all.v7.1.symbols.gmt'
    c7_path <- 'data/TMA36_project/RNA_Seq/GSEA/msigdb_v7.1_GMTs/c7.all.v7.1.symbols.gmt'
    msig_path <- 'data/TMA36_project/RNA_Seq/GSEA/msigdb_v7.1_GMTs/msigdb.v7.1.symbols.gmt'

    # For REACTOME pathways (only work with ENTREZID)

    fgsea_fixed <- function(pthw_path, reactome=FALSE, ranks=ranks){
        if(reactome){
            library(org.Hs.eg.db)
            hs <- org.Hs.eg.db
            my.symbols <- c(names(ranks))
            entrz <- AnnotationDbi::select(hs, 
                keys = my.symbols,
                columns = c("ENTREZID", "SYMBOL"),
                keytype = "SYMBOL")
            entrz <- na.omit(entrz)
            ranks <- ranks[entrz$SYMBOL]
            names(ranks) <- entrz$ENTREZID
            pth <- reactomePathways(names(ranks))
        }else{
            pth <- gmtPathways(pthw_path)
        }
        res <- fgseaMultilevel(pathways=pth, stats=ranks, eps = 0, nPermSimple = 10000) %>% #multilevel
        #res <- fgsea(pathways=pth, stats=ranks, eps = 0) %>% 
                    as_tibble %>%
                    arrange(padj, desc(abs(NES)))
        res$state <- ifelse(res$NES > 0, "up", "down")
        res$leadingEdge <- sapply(res$leadingEdge, . %>% {
          str_c(., collapse = " ")})

        res
    }



    res_hm <- fgsea_fixed(pthw_path=hallmark_path, reactome=FALSE, ranks=ranks)
    res_c1 <- fgsea_fixed(pthw_path=c1_path, reactome=FALSE, ranks=ranks)
    res_c2 <- fgsea_fixed(pthw_path=c2_path, reactome=FALSE, ranks=ranks)
    res_c3 <- fgsea_fixed(pthw_path=c3_path, reactome=FALSE, ranks=ranks)
    res_c4 <- fgsea_fixed(pthw_path=c4_path, reactome=FALSE, ranks=ranks)
    res_c5 <- fgsea_fixed(pthw_path=c5_path, reactome=FALSE, ranks=ranks)
    res_c6 <- fgsea_fixed(pthw_path=c6_path, reactome=FALSE, ranks=ranks)
    res_c7 <- fgsea_fixed(pthw_path=c7_path, reactome=FALSE, ranks=ranks)
    res_msg <- fgsea_fixed(pthw_path=msig_path, reactome=FALSE, ranks=ranks)
    res_rtm <- fgsea_fixed(reactome=TRUE, ranks=ranks)

    fgsea_res <- list(res_hm=res_hm, res_c1=res_c1, res_c2=res_c2, 
        res_c3=res_c3, res_c4=res_c4, res_c5=res_c5, res_c6=res_c6,
        res_c7=res_c7, res_msg=res_msg, res_rtm=res_rtm)

    fgsea_res
}

##########################################################################
# Pathway analysis with KEGG pathways and GO
##########################################################################
kegg_go <- function(DE_res, kegg = TRUE, GO = FALSE){
    
    library(gage)
    library(pathview)
    library(gageData)
    data(kegg.sets.hs)
    data(go.sets.hs)

    # Function to put results in a df
    ls2df <- function(pth_res){
        pth_res <- pth_res[['greater']] %>%
                as.data.frame() %>%
                na.omit()
        pth_res$state <- ifelse(pth_res$stat.mean > 0, "up", "down")
        pth_res
    }  

    # setting up the input
    res_df=DE_res$res_df
    res_df$ENSEMBL <- sapply(strsplit(res_df$gene, "\\."), "[[", 1)

    library(org.Hs.eg.db)
    hs <- org.Hs.eg.db
    ensmbl <- res_df$ENSEMBL
    entrz <- AnnotationDbi::select(hs, 
        keys = ensmbl,
        columns = c("ENTREZID", "ENSEMBL"),
        keytype = "ENSEMBL") %>% na.omit()
    entrz <- merge(res_df, entrz, by='ENSEMBL')

    fc <- entrz$log2FoldChange
    names(fc) <- entrz$ENTREZID

    # running the KEGG pathway analysis
    if(kegg){
        keggres <- gage(fc, gsets=kegg.sets.hs, same.dir=FALSE)
        keggres <- ls2df(keggres)
        return(keggres)
    }

    # Running GO 
    if(GO){
        gores = gage(fc, gsets=go.sets.hs, same.dir=FALSE)
        gores <- ls2df(gores) 
        return(gores)       
    }
    #https://www.r-bloggers.com/2015/12/tutorial-rna-seq-differential-expression-pathway-analysis-with-sailfish-deseq2-gage-and-pathview/
}

# plots
heatmap_200 <- function(res_df, vsd_mat, meta_data, pData_rnaseq, n_genes=200,
    pval_cutoff = 0.05, l2fc_cutoff=1.5, ha_custom=NULL, row_km=2, scale_mat = F,
    lowhigh = FALSE){
    res_df <- data.frame(res_df) %>%
        na.omit() %>%
        filter(abs(log2FoldChange) > l2fc_cutoff) %>%
        filter(pvalue < pval_cutoff) %>%
        arrange(pvalue)

    if(nrow(res_df)>n_genes){
        gene_list <- res_df %>%
            head(n_genes) %>%
            dplyr::select(gene) %>%
            pull()
    }else{
        gene_list <- res_df %>%
            dplyr::select(gene) %>%
            pull()        
    }

    # Use normalized data (vsd_mat) for plot
    filtered_res <- data.frame(vsd_mat) %>%
      filter(gene %in% gene_list) %>%
      as.data.frame()
    rownames(filtered_res) <- make.names(filtered_res$symbol, unique=TRUE)
    filtered_res$gene <- NULL
    filtered_res$symbol <- NULL
    
    if(lowhigh){
      ha = HeatmapAnnotation(
        Condition = as.factor(meta_data$Condition),
        col = list(Condition = c("low"="#3498DB", "high"="#EC7063")),
        simple_anno_size = unit(0.5, "cm")
      )
    }else{
      ha = HeatmapAnnotation(
        Condition = as.factor(meta_data$Condition),
        simple_anno_size = unit(0.5, "cm")
      )
    }

    if(!is.null(ha_custom)){
        ha <- ha_custom
    }

    if(scale_mat){
        filtered_res <- t(scale(t(as.matrix(filtered_res))))
    }
    print(Heatmap(filtered_res, name = "mat", 
        #column_km = 2, 
      #row_km = row_km,
      column_split =as.factor(meta_data$Condition),
      heatmap_legend_param = list(color_bar = "continuous"), 
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 8), top_annotation = ha))
}

volcano_plot <- function(res_df, gene=NULL, p_title, pCutoff=0.05, FCcutoff=1.5){

    if (!is.null(gene)){
        k <- which(res_df$gene %in% gene)
        res_df <- res_df[-k,]
    }
    res_df <- na.omit(res_df) #remove NA, these are outliers
    print(EnhancedVolcano(res_df,
        lab = rownames(res_df),
        x = 'log2FoldChange',
        y = 'pvalue',
        pCutoff=pCutoff,
        FCcutoff=FCcutoff,
        xlim = c(-5,5),
        pointSize=2,
        labSize=4,
        title = p_title))
}

fgsea_plot <- function(fgsea_res, pathways_title, cutoff = 0.05, 
    max_pathways = 30, condition_name, down_color='lightblue', up_color='#DC143C', 
    write_pathways=FALSE){

        color_levels <- function(fgsea_res) {
            colors <- c()
            if (any(fgsea_res$state == "down")) {
              colors <- c(colors, down_color)
            }
            if (any(fgsea_res$state == "up")) {
              colors <- c(colors, up_color)
            }
            colors
        }

        # Add * code for p vals
        fgsea_res$pvlabel <- '*'
        fgsea_res$pvlabel[which(fgsea_res$padj <0.01 & fgsea_res$padj>0.001)] <- '**'
        fgsea_res$pvlabel[which(fgsea_res$padj<0.001)] <- '***'
        
        # Add cols for geom_text features (hjust and y)
        fgsea_res$y <- 0.05
        fgsea_res$y[which(fgsea_res$NES>0)] <- -0.05
        fgsea_res$hjust <- 'left'
        fgsea_res$hjust[which(fgsea_res$NES>0)] <- 'right'

        if (!is.null(cutoff)) {
            fgsea_res <- fgsea_res %>% filter(padj < cutoff)
        }
        
        curated_pathways <- fgsea_res %>%
                arrange(desc(abs(NES))) %>%
                dplyr::slice(1:max_pathways)
        
        curated_pathways['leadingEdge'] <- NULL
        
        print(ggplot(curated_pathways, aes(reorder(pathway, NES), NES)) +
            geom_col(aes(fill = state), width = 0.6, color = "black") +
            scale_size_manual(values = c(0, 1), guide = "none") +
            # geom_text(
            #           aes(pathway, label = pathway, vjust=0.5, hjust=hjust, y=y),
            #           inherit.aes = TRUE, size=2)+
            geom_label(aes(label = pvlabel), size = 3, alpha = 0.75) +
            coord_flip() +
            labs(
                x = 'Pathway', 
                y = "Normalized Enrichment Score",
                title = str_c(pathways_title, " pathways: ", condition_name)
                #subtitle = str_c("(Cutoff: p.adj <", cutoff, ")")
            ) +
            theme_bw() +
            theme(axis.text = element_text(size = 10, color='black'),
                  #axis.text.y = element_blank(),
                  axis.title = element_text(size = 12, color='black'),
                  plot.title = element_blank(),#element_text(size = 15, hjust = 0.5),
                  axis.title.y = element_blank(),
                  #axis.ticks.y = element_blank(),
                  legend.position = "none")+
            scale_fill_manual(values = color_levels(curated_pathways)))

        fgsea_res <- fgsea_res %>% 
                dplyr::select(-leadingEdge, -ES, -y, -hjust) %>% 
                arrange(desc(abs(NES)))
                #DT::datatable())
        if(write_pathways) {
          write.csv(fgsea_res, file = paste0('Reports/rnaseq/tables_paper/', 
                                             'fgsea_', pathways_title, '_',
                                             gsub(' ', '_', condition_name), '.csv'))
        }

        #knitr::kable(fgsea_res)
        DT::datatable(fgsea_res, options = list(autoWidth = FALSE, scrollX=TRUE))
}


keggGO_plot <- function(keggGO_res, pathways_title, cutoff = 0.05, 
    max_pathways = 30, condition_name, pval_colnm){

        color_levels <- function(fgsea_res) {
            colors <- c()
            if (any(keggGO_res$state == "down")) {
              colors <- c(colors, "lightblue")
            }
            if (any(keggGO_res$state == "up")) {
              colors <- c(colors, "#DC143C")
            }
            colors
        }

        # Add * code for p vals
        keggGO_res$pvlabel <- '*'
        keggGO_res$pvlabel[which(keggGO_res[,pval_colnm] <0.01 & keggGO_res[,pval_colnm]>0.001)] <- '**'
        keggGO_res$pvlabel[which(keggGO_res[,pval_colnm]<0.001)] <- '***'

        # Add column for pathway and pathway_id
        pathway_id <- sapply(strsplit(rownames(keggGO_res), " "), "[[", 1)
        pathway<- sub("^(?:\\S+\\s+)", "\\1", rownames(keggGO_res), perl = TRUE)
        keggGO_res <- cbind(pathway_id, pathway, keggGO_res)
        rownames(keggGO_res)<-NULL

        if (!is.null(cutoff)) {
            keggGO_res$pvselect <- keggGO_res[,pval_colnm]
            keggGO_res <- keggGO_res %>% filter(pvselect < cutoff)
        }
        
        curated_pathways <- keggGO_res %>%
                arrange(desc(abs(stat.mean))) %>%
                dplyr::slice(1:max_pathways)

        print(ggplot(curated_pathways, aes(reorder(pathway, stat.mean), stat.mean)) +
            geom_col(aes(fill = state), width = 0.5, color = "black") +
            scale_size_manual(values = c(0, 1), guide = "none") +
            geom_label(aes(label = pvlabel), size = 3, alpha = 0.75) +
            coord_flip() +
            labs(
                x = '', 
                y = "stats.mean",
                title = str_c(pathways_title, ": ", condition_name),
                subtitle = str_c("(Cutoff: ", pval_colnm, " < ", cutoff, ")")
            ) +
            theme_bw() +
            scale_fill_manual(values = color_levels(curated_pathways)))

        keggGO_res <- keggGO_res %>% 
                dplyr::select(-exp1, -pvselect) %>% 
                arrange(desc(abs(stat.mean)))
                #DT::datatable())

        #knitr::kable(keggGO_res)
        DT::datatable(keggGO_res,options = list(autoWidth = FALSE, scrollX=TRUE))
}



hm_genes <- function(ls_preprocessed, mat_name, hla_ids, hla_gsym, scale_mat = TRUE, 
                     cond_colname = "n_op2", delete_group = TRUE, delete_group_name = 'int'){
  
  ls_preprocessed$CDE['v'] <- ls_preprocessed$CDE[cond_colname]
  mat <- as.matrix(ls_preprocessed[[mat_name]])
  if(delete_group){
    int_idx <- which(ls_preprocessed$CDE$v == delete_group_name)
    mat <- mat[,-int_idx]
    ls_preprocessed$CDE <- ls_preprocessed$CDE[-int_idx,]
  }
  
  mat <- mat[hla_ids,]
  ha = HeatmapAnnotation(
    Condition = as.factor(ls_preprocessed$CDE$v),
    SILA = ls_preprocessed$CDE$SILA,
    Stage = ls_preprocessed$CDE$Stages_simplified,
    simple_anno_size = unit(0.5, "cm")
  )
  if(scale_mat){
    mat <- t(scale(t(mat)))
  }
  rownames(mat) <- hla_gsym
  Heatmap(mat, name = "mat", 
          column_split = as.factor(ls_preprocessed$CDE$v),
          heatmap_legend_param = list(color_bar = "continuous"), 
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 8), top_annotation = ha)
  
}


pca_hla <- function(ls_preprocessed, mat_name, hla_ids, scale_mat = TRUE, 
                    cond_colname = "n_op2", delete_group = TRUE, delete_group_name = 'int'){
  
  ls_preprocessed$CDE['v'] <- ls_preprocessed$CDE[cond_colname]
  
  mat <- data.frame(t(ls_preprocessed[[mat_name]][hla_ids,]))
  mat <- cbind(mat, 'SILA'=ls_preprocessed$CDE$SILA,
               "v"=ls_preprocessed$CDE$v)
  
  if(delete_group){
    int_idx <- which(ls_preprocessed$CDE$v == delete_group_name)
    mat <- mat[-int_idx,]
  }
  
  x <- mat[,-which(colnames(mat) %in% c("SILA","v"))]
  
  if(scale_mat){
    x <- scale(x)
  }
  pca_hla <- cbind(data.frame(prcomp(x)$x), 'SILA'=mat$SILA,
                   "v"=mat$v)
  ggplot(pca_hla, aes(PC1, PC2, color = v)) +
    geom_point()
  
}


compare_hla <- function(ls_preprocessed, mat_name, hla_ids, hla_gsym, BPplot = T, 
                        cond_colname = "n_op2", delete_group = TRUE, delete_group_name = 'int'){
  
  ls_preprocessed$CDE['v'] <- ls_preprocessed$CDE[cond_colname]
  
  if(length(hla_ids)<2){
    mat <- data.frame(ls_preprocessed[[mat_name]][hla_ids,])
    colnames(mat) <- hla_ids
  }else{
    mat <- data.frame(t(ls_preprocessed[[mat_name]][hla_ids,]))
  }
  
  mat <- cbind(mat, 'SILA'=ls_preprocessed$CDE$SILA,
               "v"=ls_preprocessed$CDE$v)
  if(delete_group){
    int_idx <- which(ls_preprocessed$CDE$v == delete_group_name)
    mat <- mat[-int_idx,]
  }
  
  sum_df <- data.frame(matrix(nrow = 1, ncol = 4))
  groups_c <- unique(mat$v)
  g1_idx <- which(mat$v == groups_c[1])
  for (i in colnames(mat)[grep("ENS",colnames(mat))]){
    g1_med <- median(mat[g1_idx,i])
    g2_med <- median(mat[-g1_idx,i])
    pval <- wilcox.test(mat[g1_idx,i], mat[-g1_idx,i])
    sum_df <- rbind(sum_df, c(i,g1_med, g2_med, pval$p.value))
  }
  sum_df <- sum_df[-1,]
  colnames(sum_df) <- c('Gene', paste0('Median.', groups_c[1]), paste0('Median.', groups_c[2]), 'p.value')
  sum_df['p.adjusted'] <- p.adjust(sum_df$p.value, method = "BH")
  sum_df <- cbind('Symbol'=hla_gsym, sum_df)
  sum_df <- sum_df %>%
    arrange(p.value) %>%
    slice(which(p.value < 0.09))
  
  # boxplots
  if(BPplot){
    x <- data.frame(mat[,sum_df$Gene])
    colnames(x) <- sum_df$Symbol
    x['v'] <- mat$v
    x <- reshape2::melt(x)
    lb <- x$v
    #lb[which(x[,'v'] == 'agg')] <- 'Aggressive'
    #lb[-which(x[,'v'] == 'agg')] <- 'Indolent'
    #levels = c('Indolent', 'Aggressive')
    x['v'] <- factor(lb, levels = groups_c)
    
    #comp =  list(c("Indolent", "Aggressive"))
    comp =  list(groups_c)
    colr = c("#3498DB", "#EC7063")
    
    if(nrow(sum_df)>1){
      ggplot(x, aes(x=v, y=value, color = v)) +
        geom_boxplot() +
        ggsignif::geom_signif(comparisons = comp, 
                              map_signif_level=TRUE) +
        facet_wrap(~variable, scales='free') +
        theme(plot.title = element_text(hjust = 0.5, size=22))+
        scale_color_manual(values=colr, name = "Behavior")       
    }else{
      #ggplot(x, aes(x=v, y=value, color = v)) +
    }
    
  }
  #knitr::kable(sum_df)
  DT::datatable(sum_df, options = list(autoWidth = FALSE, scrollX=TRUE))
  
}


survival_plot <- function(CDE, group_colname, group_levels, delete_group=FALSE, delete_group_name,
                          survival_type = 'OS', legend_labs, legend_title, plot_rmst=FALSE,
                          rmst2_tau=4, col_pal= c("#3498DB", "#EC7063")){
  # generate survival data
  surv_data <- CDE %>%
    mutate(Recurrence_Date = ifelse(is.na(Recurrence_Date), LDKA, Recurrence_Date),
           Progression_Date = ifelse(is.na(Progression_Date), LDKA, Progression_Date),
           OS_yrs = lubridate::time_length(difftime(as.Date(LDKA), as.Date(Thoracotomy_Date)), "years"), #Diagnosis_Date
           RFS_yrs = lubridate::time_length(difftime(as.Date(Recurrence_Date), as.Date(Thoracotomy_Date)), "years"),
           PFS_yrs = lubridate::time_length(difftime(as.Date(Progression_Date), as.Date(Thoracotomy_Date)), "years"),
           Death_st = ifelse(Death_st=='Yes', 1, 0),
           Recurrence_st = ifelse(Recurrence_st=='Yes', 1, 0),
           Progression_st = ifelse(Progression_st=='Yes', 1, 0),
           DRP_st = ifelse(DRP_st=='Yes', 1, 0))%>%
    mutate(., DRP_yrs=apply(.[,c('OS_yrs', 'RFS_yrs', 'PFS_yrs')], 1, min))%>%
    dplyr::select(., pt_ID, !!group_colname, SILA, OS_yrs, RFS_yrs, PFS_yrs, DRP_yrs, Death_st, Recurrence_st, Progression_st, DRP_st)
  
  surv_data['group'] <- factor(surv_data[,group_colname], levels = group_levels)
  
  if(length(group_levels)>2){
    col_pal <- RColorBrewer::brewer.pal(length(group_levels), 'Dark2')
    dt <- surv_data[-which(surv_data$group==group_levels[2]),]
    #col_pal <- c("#3498DB", "#888a87", "#EC7063")
  }else{
    dt <- surv_data
    col_pal <- col_pal
  }
  
  time = dt$OS_yrs
  status = dt$Death_st
  
  if(delete_group){
    surv_data <- surv_data[-which(surv_data$group == delete_group_name)]
  }
  
  # create survival object
  if(survival_type == 'RFS'){
    km_Group_fit <- survfit(Surv(RFS_yrs, Recurrence_st) ~ group, data = surv_data)
    title_plt <- paste('Recurrence Free Survival by', legend_title)
    time = dt$RFS_yrs
    status = dt$Recurrence_st
  }else if(survival_type == 'PFS'){
    km_Group_fit <- survfit(Surv(PFS_yrs, Progression_st) ~ group, data = surv_data)
    title_plt <- paste('Progression Free Survival by', legend_title)
    time = dt$PFS_yrs
    status = dt$Progression_st
  }else if(survival_type == 'DRP'){
    km_Group_fit <- survfit(Surv(DRP_yrs, DRP_st) ~ group, data = surv_data)
    title_plt <- paste('Overall/Recurrence/Progression Free Survival by', legend_title)
    time = dt$DRP_yrs
    status = dt$DRP_st
  }else{
    km_Group_fit <- survfit(Surv(OS_yrs, Death_st) ~ group, data = surv_data)
    title_plt <- paste('Overall Survival by', legend_title)
    time = dt$OS_yrs
    status = dt$Death_st
  }
  
  # Plot survival curves
  p <- ggsurvplot(km_Group_fit, data=surv_data, palette = col_pal,
                  fun = 'pct', ggtheme = theme_bw(), risk.table = TRUE, 
                  tables.theme = theme_cleantable(), surv.median.line = "hv", 
                  xlab = 'Years', ylab = 'Survival Probability (%)',
                  title = title_plt, pval = TRUE, 
                  legend.title = legend_title, legend.labs = legend_labs)
  
  print(p)
  #print(p$plot)
  #print(p$table)
  survival_ls <- list(surv_data = surv_data, 
                      km_Group_fit = km_Group_fit)
  # RMST calculations 
  if(plot_rmst){
    arm = ifelse(dt$group==group_levels[1], 1,0)
    rmst_calc = rmst2(time, status, arm, tau=rmst2_tau)
    plot(rmst_calc, xlab="Years", ylab="Probability", density=60)
    survival_ls$rmst_calc <- rmst_calc
  }
  
  
  return(survival_ls)
  
}
