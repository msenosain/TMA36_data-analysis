# Import source code
source("src/cytof/20_ClustAnnot_functions.R")

#########################################################
# Data without CD90
#########################################################

# Load data
load("data/TMA36_project/CyTOF/processed/Data_paper2/woCD90/woCD90.RData")

# Clustering
dt_cl <- clustering(big_df, n_clusters = 10, iterations = 200, seed = 45) 
# 15 20 29 31 37 40

# Clustering Evaluation
hm_data <- ClusterEval_data(dt_cl, eval_type = 'heatmap')
ClusterEval_plot(hm_data, data_type = 'heatmap')

tsne_data <- ClusterEval_data(dt_cl, eval_type = 'tSNE', sample_size = 25000)
ClusterEval_plot(tsne_data, data_type = 'tSNE')

umap_data <- ClusterEval_data(dt_cl, eval_type = 'UMAP', sample_size = 25000)
ClusterEval_plot(umap_data, data_type = 'UMAP')

table(dt_cl$cluster)/nrow(dt_cl)*100

# Cluster annotation
epithelial <- c(6,7,9)
endothelial <- c(8)
mesenchymal <- c(4)
immune <- c(1,2,3,10) # for manuscript 1 cluster 10 is epithelial
nothing <- c(5)
#unknown <- c(6)

ct_ls <- list('Epithelial'= epithelial, 'Endothelial'=endothelial, 
                   'Mesenchymal'=mesenchymal, 'Immune'=immune, 'Nothing'=nothing)

annot_df_wo <- ClusterAnnotation(data = big_df, df_cluster = dt_cl, 
    ls_annotation = ct_ls, annotation_col = 'cell_type', cl_delete = T, 
    cl_delete_name = 'Nothing')

ref_wo <- ref

save(hm_data, tsne_data, umap_data, annot_df_wo, ref_wo, file = 'data/TMA36_project/CyTOF/processed/Data_paper2/woCD90/majorcelltypes.RData')


# Do this for both, wCD90 woCD90 separately and then merge.
# Subset cell types for further subtyping

#########################################################
# Data with CD90
#########################################################

# Load data
load("data/TMA36_project/CyTOF/processed/Data_paper2/withCD90/withCD90.RData")

# Clustering
dt_cl <- clustering(big_df, n_clusters = 12, iterations = 200, seed = 45)
# 15 20 29 31 37 40 49

# Clustering Evaluation
hm_data <- ClusterEval_data(dt_cl, eval_type = 'heatmap')
ClusterEval_plot(hm_data, data_type = 'heatmap')

tsne_data <- ClusterEval_data(dt_cl, eval_type = 'tSNE', sample_size = 25000)
ClusterEval_plot(tsne_data, data_type = 'tSNE')

umap_data <- ClusterEval_data(dt_cl, eval_type = 'UMAP', sample_size = 25000)
ClusterEval_plot(umap_data, data_type = 'UMAP')

table(dt_cl$cluster)/nrow(dt_cl)*100

# Cluster annotation
epithelial <- c(1,2)
endothelial <- c(7)
fibroblasts <- c(9)
mesenchymal <- c(11)
immune <- c(3,4,5,6,8,12)
nothing <- c(10)

ct_ls <- list('Epithelial'= epithelial, 'Endothelial'=endothelial, 
                   'Fibroblasts' = fibroblasts, 'Mesenchymal'=mesenchymal, 
                   'Immune'=immune, 'Nothing'=nothing)

annot_df <- ClusterAnnotation(data = big_df, df_cluster = dt_cl, 
    ls_annotation = ct_ls, annotation_col = 'cell_type', cl_delete = T, 
    cl_delete_name = 'Nothing')

save(hm_data, tsne_data, umap_data, annot_df, ref, file = 'data/TMA36_project/CyTOF/processed/Data_paper2/withCD90/majorcelltypes.RData')


#########################################################
# Combining both w and wo CD90
#########################################################

annot_df <- change_colname(annot_df, 'cell_type', 'cell_type_A')
annot_df_wo <- change_colname(annot_df_wo, 'cell_type', 'cell_type_A')

annot_df['cell_type_B'] <- annot_df$cell_type_A
annot_df_wo['cell_type_B'] <- annot_df_wo$cell_type_A

annot_df <- edit_names(annot_df, col_name = 'cell_type_B', var_oldname = 'Fibroblasts',
                var_newname = 'Fib_Mesenchymal')
annot_df <- edit_names(annot_df, col_name = 'cell_type_B', var_oldname = 'Mesenchymal',
                var_newname = 'Fib_Mesenchymal')

annot_df_wo <- edit_names(annot_df_wo, col_name = 'cell_type_B', var_oldname = 'Mesenchymal',
                var_newname = 'Fib_Mesenchymal')

annot_df <- rbind(annot_df_wo, annot_df)
ref <- rbind(ref_wo, ref)

dir.create("data/TMA36_project/CyTOF/processed/Data_paper2/both")

save(annot_df, ref, file = "data/TMA36_project/CyTOF/processed/Data_paper2/both/majorcelltypes_merged.RData")


# Immune cells
##########################################
load("data/TMA36_project/CyTOF/processed/Data_paper2/both/majorcelltypes_merged.RData")

# Clustering for Immune cells
ct_immune <- subset(annot_df, cell_type_A == 'Immune')

DetermineNumberOfClusters(ct_immune, k_max=15, plot=T,smooth=0.2,
                                      iter.max=200, seed = 45)

imm_cl <- clustering(ct_immune, n_clusters = 8, iterations = 200, seed = 55)
# 30 34 44 46 47

# Clustering Evaluation
hm_data <- ClusterEval_data(imm_cl, eval_type = 'heatmap')
ClusterEval_plot(hm_data, data_type = 'heatmap')
table(imm_cl$cluster)/nrow(imm_cl)*100

umap_data <- ClusterEval_data(imm_cl, eval_type = 'UMAP', sample_size = 25000)
ClusterEval_plot(umap_data, data_type = 'UMAP')

# Cluster annotation
myeloid <- c(3,7) # 3 7
th_cells <- c(6) # 6
tc_cells <- c(2) # 2
dnt_cells <- c(4,5) # 4 5
other_immune <- c(1,8) # 1 8

im_ls <- list('Myeloid'= myeloid, 'Th_cells'=th_cells, 'Tc_cells'=tc_cells, 
    'DNT_cells'= dnt_cells, 'Other_immune'=other_immune)

annot_imm <- ClusterAnnotation(data = ct_immune, df_cluster = imm_cl, 
    ls_annotation = im_ls, annotation_col = 'subtype')

dir <- 'data/TMA36_project/CyTOF/processed/Data_paper2/both'

save(hm_data, umap_data, annot_imm, file = file.path(dir,'immunesubtypes.RData'))



# Other cell types
##########################################
load("data/TMA36_project/CyTOF/processed/Data_paper2/both/immunesubtypes.RData")
load("data/TMA36_project/CyTOF/processed/Data_paper2/both/majorcelltypes_merged.RData")

annot_imm <- change_colname(annot_imm, 'subtype', 'subtype_A')

annot_imm['subtype_B'] <- annot_imm$subtype_A

ct_epi <- subset(annot_df, cell_type_A == 'Epithelial')
ct_epi[,'subtype_A'] <- ct_epi$cell_type_A
ct_epi[,'subtype_B'] <- ct_epi$cell_type_A

ct_mes <- subset(annot_df, cell_type_A == 'Mesenchymal')
ct_mes[,'subtype_A'] <- ct_mes$cell_type_A

ct_fib <- subset(annot_df, cell_type_A == 'Fibroblasts')
ct_fib[,'subtype_A'] <- ct_fib$cell_type_A

ct_fibm <- rbind(ct_mes, ct_fib)
ct_fibm[,'subtype_B'] <- ct_fibm$cell_type_B

ct_endo <- subset(annot_df, cell_type_A == 'Endothelial')
ct_endo[,'subtype_A'] <- ct_endo$cell_type_A
ct_endo[,'subtype_B'] <- ct_endo$cell_type_A


annot_df <- rbind(annot_imm, ct_epi, ct_fibm, ct_endo)

dir <- 'data/TMA36_project/CyTOF/processed/Data_paper2/both'

save(ref, annot_df, file = file.path(dir,'cellsubtypes.RData'))


##########################################
# Clustering cell types
##########################################
load("data/TMA36_project/CyTOF/processed/Data_paper2/both/cellsubtypes.RData")

# Clustering Epithelial cells
##########################################
ct_epi <- subset(annot_df, cell_type_A == 'Epithelial')

DetermineNumberOfClusters(ct_epi, k_max=15, plot=T,smooth=0.2,
                                      iter.max=200, seed = 45)

epi_cl <- clustering(ct_epi, n_clusters = 6, iterations = 200, seed = 42)
# 15 18 19 21 22 23 26 28 31 33 35 37 38 40 41 50 51

# Clustering Evaluation
hm_data <- ClusterEval_data(epi_cl, eval_type = 'heatmap')
ClusterEval_plot(hm_data, data_type = 'heatmap')

table(epi_cl$cluster)/nrow(epi_cl)*100

umap_data <- ClusterEval_data(epi_cl, eval_type = 'UMAP', sample_size = 25000)
ClusterEval_plot(umap_data, data_type = 'UMAP')

# Cluster annotation
epi_ls <- list()
for(i in 1:6){
    epi_ls[i] <- c(i)
}
names(epi_ls) <- paste0('ECC_',rep(1:6))

annot_epi <- ClusterAnnotation(data = ct_epi, df_cluster = epi_cl, 
    ls_annotation = epi_ls, annotation_col = 'clusters_A')
annot_epi['clusters_B'] <- annot_epi$clusters_A

dir <- 'data/TMA36_project/CyTOF/processed/Data_paper2/both'
save(hm_data, umap_data, annot_epi, file = file.path(dir,'epiclusters.RData'))


# Clustering CD4 T cells
##########################################
# 17 18 19 20 21 22 24 25 26 28 31 33 35 42 45 48 50 # immune
load("data/TMA36_project/CyTOF/processed/Data_paper2/both/cellsubtypes.RData")

ct_cd4 <- subset(annot_df, subtype_A == 'Th_cells')

DetermineNumberOfClusters(ct_cd4, k_max=10, plot=T,smooth=0.2,
                                      iter.max=200, seed = 45)

cd4_cl <- clustering(ct_cd4, n_clusters = 3, iterations = 200, seed = 50)
# 18 19 20 21 22 26 28 31 35 42 50

# Clustering Evaluation
hm_data <- ClusterEval_data(cd4_cl, eval_type = 'heatmap')
ClusterEval_plot(hm_data, data_type = 'heatmap')

table(cd4_cl$cluster)/nrow(cd4_cl)*100

umap_data <- ClusterEval_data(cd4_cl, eval_type = 'UMAP', sample_size = 10000)
ClusterEval_plot(umap_data, data_type = 'UMAP')

# Cluster annotation
cd4_ls <- list()
for(i in 1:3){
    cd4_ls[i] <- c(i)
}
names(cd4_ls) <- paste0('Th_',rep(1:3))

annot_cd4 <- ClusterAnnotation(data = ct_cd4, df_cluster = cd4_cl, 
    ls_annotation = cd4_ls, annotation_col = 'clusters_A')
annot_cd4['clusters_B'] <- annot_cd4$clusters_A

dir <- 'data/TMA36_project/CyTOF/processed/Data_paper2/both'
save(hm_data, umap_data, annot_cd4, file = file.path(dir,'cd4clusters.RData'))



# Clustering CD8 T cells
##########################################
ct_cd8 <- subset(annot_df, subtype_A == 'Tc_cells')

DetermineNumberOfClusters(ct_cd8, k_max=10, plot=T,smooth=0.2,
                                      iter.max=200, seed = 45)

cd8_cl <- clustering(ct_cd8, n_clusters = 3, iterations = 200, seed = 50)
# 18 19 20 21 22 26 28 31 35 42 50

# Clustering Evaluation
hm_data <- ClusterEval_data(cd8_cl, eval_type = 'heatmap')
ClusterEval_plot(hm_data, data_type = 'heatmap')

table(cd8_cl$cluster)/nrow(cd8_cl)*100

umap_data <- ClusterEval_data(cd8_cl, eval_type = 'UMAP', sample_size = 10000)
ClusterEval_plot(umap_data, data_type = 'UMAP')

# Cluster annotation
cd8_ls <- list()
for(i in 1:3){
    cd8_ls[i] <- c(i)
}
names(cd8_ls) <- paste0('Tc_',rep(1:3))

annot_cd8 <- ClusterAnnotation(data = ct_cd8, df_cluster = cd8_cl, 
    ls_annotation = cd8_ls, annotation_col = 'clusters_A')
annot_cd8['clusters_B'] <- annot_cd8$clusters_A

dir <- 'data/TMA36_project/CyTOF/processed/Data_paper2/both'
save(hm_data, umap_data, annot_cd8, file = file.path(dir,'cd8clusters.RData'))

# Clustering DN T cells
##########################################
ct_dnt <- subset(annot_df, subtype_A == 'DNT_cells')

DetermineNumberOfClusters(ct_dnt, k_max=10, plot=T,smooth=0.2,
                                      iter.max=200, seed = 45)

dnt_cl <- clustering(ct_dnt, n_clusters = 3, iterations = 200, seed = 50)
# 18 19 20 21 22 26 28 31 35 42 50

# Clustering Evaluation
hm_data <- ClusterEval_data(dnt_cl, eval_type = 'heatmap')
ClusterEval_plot(hm_data, data_type = 'heatmap')

table(dnt_cl$cluster)/nrow(dnt_cl)*100

umap_data <- ClusterEval_data(dnt_cl, eval_type = 'UMAP', sample_size = 10000)
ClusterEval_plot(umap_data, data_type = 'UMAP')

# Cluster annotation
dnt_ls <- list()
for(i in 1:3){
    dnt_ls[i] <- c(i)
}
names(dnt_ls) <- paste0('DNT_',rep(1:3))

annot_dnt <- ClusterAnnotation(data = ct_dnt, df_cluster = dnt_cl, 
    ls_annotation = dnt_ls, annotation_col = 'clusters_A')
annot_dnt['clusters_B'] <- annot_dnt$clusters_A

dir <- 'data/TMA36_project/CyTOF/processed/Data_paper2/both'
save(hm_data, umap_data, annot_dnt, file = file.path(dir,'dntclusters.RData'))

# Clustering Myeloid cells
##########################################
ct_mye <- subset(annot_df, subtype_A == 'Myeloid')

DetermineNumberOfClusters(ct_mye, k_max=10, plot=T,smooth=0.2,
                                      iter.max=200, seed = 45)

mye_cl <- clustering(ct_mye, n_clusters = 3, iterations = 200, seed = 50)
# 18 19 20 21 22 26 28 31 35 42 50

# Clustering Evaluation
hm_data <- ClusterEval_data(mye_cl, eval_type = 'heatmap')
ClusterEval_plot(hm_data, data_type = 'heatmap')

table(mye_cl$cluster)/nrow(mye_cl)*100

umap_data <- ClusterEval_data(mye_cl, eval_type = 'UMAP', sample_size = 10000)
ClusterEval_plot(umap_data, data_type = 'UMAP')

# Cluster annotation
mye_ls <- list()
for(i in 1:3){
    mye_ls[i] <- c(i)
}
names(mye_ls) <- paste0('Myeloid_',rep(1:3))

annot_mye <- ClusterAnnotation(data = ct_mye, df_cluster = mye_cl, 
    ls_annotation = mye_ls, annotation_col = 'clusters_A')
annot_mye['clusters_B'] <- annot_mye$clusters_A

dir <- 'data/TMA36_project/CyTOF/processed/Data_paper2/both'
save(hm_data, umap_data, annot_mye, file = file.path(dir,'myeclusters.RData'))

# Clustering Other Immune cells
##########################################
ct_oi <- subset(annot_df, subtype_A == 'Other_immune')

DetermineNumberOfClusters(ct_oi, k_max=10, plot=T,smooth=0.2,
                                      iter.max=200, seed = 45)

oi_cl <- clustering(ct_oi, n_clusters = 4, iterations = 200, seed = 50)
# 18 19 20 21 22 26 28 31 35 42 50

# Clustering Evaluation
hm_data <- ClusterEval_data(oi_cl, eval_type = 'heatmap')
ClusterEval_plot(hm_data, data_type = 'heatmap')

table(oi_cl$cluster)/nrow(oi_cl)*100

umap_data <- ClusterEval_data(oi_cl, eval_type = 'UMAP', sample_size = 10000)
ClusterEval_plot(umap_data, data_type = 'UMAP')

# Cluster annotation
oi_ls <- list()
for(i in 1:4){
    oi_ls[i] <- c(i)
}
names(oi_ls) <- paste0('OtherI_',rep(1:4))

annot_oi <- ClusterAnnotation(data = ct_oi, df_cluster = oi_cl, 
    ls_annotation = oi_ls, annotation_col = 'clusters_A')
annot_oi['clusters_B'] <- annot_oi$clusters_A

dir <- 'data/TMA36_project/CyTOF/processed/Data_paper2/both'
save(hm_data, umap_data, annot_oi, file = file.path(dir,'otherimmclusters.RData'))

# Clustering Fib_Mesenchymal cells
##########################################
load("data/TMA36_project/CyTOF/processed/Data_paper2/both/cellsubtypes.RData")

ct_fmes <- subset(annot_df, cell_type_B == 'Fib_Mesenchymal')

DetermineNumberOfClusters(ct_fmes, k_max=10, plot=T,smooth=0.2,
                                      iter.max=200, seed = 45)

fmes_cl <- clustering(ct_fmes, n_clusters = 5, iterations = 200, seed = 50)
# 18 19 21 22 26 27 28 33 35 38 41 42 43 48 50 51

# Clustering Evaluation
hm_data <- ClusterEval_data(fmes_cl, eval_type = 'heatmap')
ClusterEval_plot(hm_data, data_type = 'heatmap')

table(fmes_cl$cluster)/nrow(fmes_cl)*100

umap_data <- ClusterEval_data(fmes_cl, eval_type = 'UMAP', sample_size = 10000)
ClusterEval_plot(umap_data, data_type = 'UMAP')

# Cluster annotation
fmes_ls <- list()
for(i in 1:5){
    fmes_ls[i] <- c(i)
}
names(fmes_ls) <- paste0('fmes_',rep(1:5))

annot_fmes <- ClusterAnnotation(data = ct_fmes, df_cluster = fmes_cl, 
    ls_annotation = fmes_ls, annotation_col = 'clusters_B')
#annot_fmes['clusters_B'] <- annot_fmes$clusters_A

####
# Clustering Mesenchymal cells
ct_mes <- subset(annot_fmes, cell_type_A == 'Mesenchymal')

DetermineNumberOfClusters(ct_mes, k_max=10, plot=T,smooth=0.2,
                                      iter.max=200, seed = 45)

mes_cl <- clustering(ct_mes, n_clusters = 4, iterations = 200, seed = 50)
# 18 19 21 22 26 27 28 33 35 38 41 42 43 48 50 51

# Clustering Evaluation
hm_data_mes <- ClusterEval_data(mes_cl, eval_type = 'heatmap')
ClusterEval_plot(hm_data_mes, data_type = 'heatmap')

table(mes_cl$cluster)/nrow(mes_cl)*100

umap_data_mes <- ClusterEval_data(mes_cl, eval_type = 'UMAP', sample_size = 10000)
ClusterEval_plot(umap_data_mes, data_type = 'UMAP')

# Cluster annotation
mes_ls <- list()
for(i in 1:4){
    mes_ls[i] <- c(i)
}
names(mes_ls) <- paste0('mes_',rep(1:4))

annot_mes <- ClusterAnnotation(data = ct_mes, df_cluster = mes_cl, 
    ls_annotation = mes_ls, annotation_col = 'clusters_A')

####
# Clustering Fibroblasts
ct_fib <- subset(annot_fmes, cell_type_A == 'Fibroblasts')

DetermineNumberOfClusters(ct_fib, k_max=10, plot=T,smooth=0.2,
                                      iter.max=200, seed = 45)

fib_cl <- clustering(ct_fib, n_clusters = 3, iterations = 200, seed = 50)
# 18 19 21 22 26 27 28 33 35 38 41 42 43 48 50 51

# Clustering Evaluation
hm_data_fib <- ClusterEval_data(fib_cl, eval_type = 'heatmap')
ClusterEval_plot(hm_data_fib, data_type = 'heatmap')

table(fib_cl$cluster)/nrow(fib_cl)*100

umap_data_fib <- ClusterEval_data(fib_cl, eval_type = 'UMAP', sample_size = 10000)
ClusterEval_plot(umap_data_fib, data_type = 'UMAP')

# Cluster annotation
fib_ls <- list()
for(i in 1:3){
    fib_ls[i] <- c(i)
}
names(fib_ls) <- paste0('fib_',rep(1:3))

annot_fib <- ClusterAnnotation(data = ct_fib, df_cluster = fib_cl, 
    ls_annotation = fib_ls, annotation_col = 'clusters_A')


annot_fmes <- rbind(annot_mes, annot_fib)
dir <- 'data/TMA36_project/CyTOF/processed/Data_paper2/both'
save(hm_data, hm_data_mes, hm_data_fib, umap_data, umap_data_mes, umap_data_fib,  annot_fmes, file = file.path(dir,'fmesclusters.RData'))

# Clustering Endothelial cells
##########################################
ct_endo <- subset(annot_df, cell_type_A == 'Endothelial')

DetermineNumberOfClusters(ct_endo, k_max=10, plot=T,smooth=0.2,
                                      iter.max=200, seed = 45)

endo_cl <- clustering(ct_endo, n_clusters = 3, iterations = 200, seed = 50)
# 18 19 21 22 26 28 31 33 35 38 41 42 50 51

# Clustering Evaluation
hm_data <- ClusterEval_data(endo_cl, eval_type = 'heatmap')
ClusterEval_plot(hm_data, data_type = 'heatmap')

table(endo_cl$cluster)/nrow(endo_cl)*100

umap_data <- ClusterEval_data(endo_cl, eval_type = 'UMAP', sample_size = 10000)
ClusterEval_plot(umap_data, data_type = 'UMAP')

# Cluster annotation
endo_ls <- list()
for(i in 1:3){
    endo_ls[i] <- c(i)
}
names(endo_ls) <- paste0('endo_',rep(1:3))

annot_endo <- ClusterAnnotation(data = ct_endo, df_cluster = endo_cl, 
    ls_annotation = endo_ls, annotation_col = 'clusters_A')
annot_endo['clusters_B'] <- annot_endo$clusters_A

dir <- 'data/TMA36_project/CyTOF/processed/Data_paper2/both'
save(hm_data, umap_data, annot_endo, file = file.path(dir,'endoclusters.RData'))

###############################################################################
# RBIND ALL CLUSTERS
load("data/TMA36_project/CyTOF/processed/Data_paper2/both/majorcelltypes_merged.RData")
remove(annot_df)
load("data/TMA36_project/CyTOF/processed/Data_paper2/both/endoclusters.RData")
load("data/TMA36_project/CyTOF/processed/Data_paper2/both/fmesclusters.RData")
load("data/TMA36_project/CyTOF/processed/Data_paper2/both/myeclusters.RData")
load("data/TMA36_project/CyTOF/processed/Data_paper2/both/dntclusters.RData")
load("data/TMA36_project/CyTOF/processed/Data_paper2/both/cd8clusters.RData")
load("data/TMA36_project/CyTOF/processed/Data_paper2/both/cd4clusters.RData")
load("data/TMA36_project/CyTOF/processed/Data_paper2/both/epiclusters.RData")
load("data/TMA36_project/CyTOF/processed/Data_paper2/both/otherimmclusters.RData")

annot_df <- rbind(annot_epi, annot_fmes, annot_endo, annot_cd4, annot_cd8, annot_dnt, annot_mye, annot_oi)

dir <- 'data/TMA36_project/CyTOF/processed/Data_paper2/both'
save(annot_df, ref, file = file.path(dir,'cellclusters.RData'))

###############################################################################
# GET Cell type percentages tables
###############################################################################
load("data/TMA36_project/CyTOF/processed/Data_paper2/both/cellclusters.RData")

# By major cell types
prcnt_celltypes <- ClassAbundanceByPt(data=annot_df, ptID_col = 'pt_ID', 
    class_col = 'cell_type_B')

# By cell subtypes
prcnt_subtypes <- ClassAbundanceByPt(data=annot_df, ptID_col = 'pt_ID', 
    class_col = 'subtype_B')

# By cell clusters
prcnt_clusters <- ClassAbundanceByPt(data=annot_df, ptID_col = 'pt_ID', 
    class_col = 'clusters_B')

dir <- 'data/TMA36_project/CyTOF/processed/Data_paper2/both'
save(prcnt_celltypes, prcnt_subtypes, prcnt_clusters, file = file.path(dir,'percent_pt.RData'))

###############################################################################
# Get protein co-expression (correlation)
###############################################################################
load("data/TMA36_project/CyTOF/processed/Data_paper2/both/cellclusters.RData")

corr_f <- function(data, rcorr_type = 'spearman', p.adjust_method = 'BH'){

    res <- Hmisc::rcorr(as.matrix(data), type = rcorr_type) #for corr plot
    # corrplot
    corrected_pvals <- p.adjust(res$P, method = p.adjust_method)
    corrected_pvals <- matrix(corrected_pvals, nrow = ncol(res$P), 
        ncol = ncol(res$P))
    colnames(corrected_pvals)<- colnames(res$P)
    rownames(corrected_pvals)<- rownames(res$P)

    res$P <- corrected_pvals

    res    
}

idx <- c(15, 17:31, 33:35, 37:51)


# Bulk
co_bulk <- corr_f(annot_df[,idx])

# By cell subtype

## Endothelial
k <- which(annot_df$cell_type_A == 'Endothelial')
co_endo <- corr_f(annot_df[k,idx])

## Epithelial
k <- which(annot_df$cell_type_A == 'Epithelial')
co_epi <- corr_f(annot_df[k,idx])

## Fib_mes
k <- which(annot_df$cell_type_B == 'Fib_Mesenchymal')
co_fmes <- corr_f(annot_df[k,idx])

## CD4 T
k <- which(annot_df$subtype_A == 'Th_cells')
co_th <- corr_f(annot_df[k,idx])

## CD8 T
k <- which(annot_df$subtype_A == 'Tc_cells')
co_tc <- corr_f(annot_df[k,idx])

## DN T
k <- which(annot_df$subtype_A == 'DNT_cells')
co_dnt <- corr_f(annot_df[k,idx])

## Myeloid
k <- which(annot_df$subtype_A == 'Myeloid')
co_mye <- corr_f(annot_df[k,idx])

## Other immune
k <- which(annot_df$subtype_A == 'Other_immune')
co_oimm <- corr_f(annot_df[k,idx])

dir <- 'data/TMA36_project/CyTOF/processed/Data_paper2/both'
save(co_bulk, co_epi, co_endo, co_fmes, co_th, co_tc, co_dnt, co_mye, co_oimm, file = file.path(dir,'protein_correlations.RData'))


#####################################################
# Controls clustering
#####################################################

load("data/TMA36_project/CyTOF/processed/Data_paper2/controls/CyTOF_ADC_controls.RData")

# Clustering
dt_cl <- clustering(big_df, n_clusters = 2, iterations = 200, seed = 45) 
# 15 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 33 34 35 37 38 39 40 41 42 43 44 45 46 47 48 50 51

# Clustering Evaluation
hm_data <- ClusterEval_data(dt_cl, eval_type = 'heatmap')
ClusterEval_plot(hm_data, data_type = 'heatmap')

tsne_data <- ClusterEval_data(dt_cl, eval_type = 'tSNE', sample_size = 25000)
ClusterEval_plot(tsne_data, data_type = 'tSNE')

umap_data <- ClusterEval_data(dt_cl, eval_type = 'UMAP', sample_size = 25000)
ClusterEval_plot(umap_data, data_type = 'UMAP')

table(dt_cl$cluster)/nrow(dt_cl)*100

# Cluster annotation
A549 <- c(2)
Ramos <- c(1)

ct_ls <- list('A549'=A549, 'Ramos'=Ramos)

annot_df_ctl <- ClusterAnnotation(data = big_df, df_cluster = dt_cl, 
    ls_annotation = ct_ls, annotation_col = 'cell_line')

dir <- 'data/TMA36_project/CyTOF/processed/Data_paper2/controls'
save(hm_data, tsne_data, umap_data, annot_df_ctl, ref, file = file.path(dir,'annotated_controls.RData'))

