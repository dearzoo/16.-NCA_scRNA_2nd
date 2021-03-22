library(Seurat)

# 0. Data ready ----
# Call pre-processed split object 
# Automated Normalize, CellCycleScoring, and SCTransform (with regressing out mt content) done
split_seurat <- readRDS(file = "objects/split_seurat.rds")


# Check which assays are stored in objects (All need to show both RNA and SCT)
for(i in 1:length(split_seurat)){
  print(names(split_seurat[[i]]@assays))
}

# index of nuc_seq samples
View(names(split_seurat))    # 5 and 6


# 1. Integration ----
# Integratino of split_seurat[[5]] and split_seurat[[6]]
# make a new list only with nuc_seq samples

split_seurat_nuc <- split_seurat[5:6]
names(split_seurat_nuc)

# Select the most variable features to use for integration (default: nfeatures = 2000)
integ_features_nuc <- SelectIntegrationFeatures(object.list = split_seurat_nuc, nfeatures = 3000)

# Prepare the SCT list object for integration
split_seurat_nuc <- PrepSCTIntegration(object.list = split_seurat_nuc,
                                   anchor.features = integ_features_nuc)

# Find best buddies - can take a while to run
integ_anchors_nuc <- FindIntegrationAnchors(object.list = split_seurat_nuc,
                                        normalization.method = "SCT",
                                        anchor.features = integ_features_nuc)

# Integrate across conditions
seurat_integrated_nuc <- IntegrateData(anchorset = integ_anchors_nuc,
                                       normalization.method = "SCT")

# save
saveRDS(seurat_integrated_nuc, file = "objects/seurat_integrated_nuc.rds")
# seurat_integrated_nuc <- readRDS(file = "objects/seurat_integrated_nuc.rds")


# 2. Reduction ----
# Run PCA
seurat_integrated_nuc <- RunPCA(object = seurat_integrated_nuc)


# re-leveling of sample column
seurat_integrated_nuc$sample <- factor(seurat_integrated_nuc$sample, levels = c("w3_47f_nuc", "w3_45f_nuc"))
table(seurat_integrated_nuc$sample)




# Plot PCA
library(ggplot2)

jpeg(filename = "figures/reduction/pca_nuc.jpeg", 
     width = 800, height = 800,
     quality = 100,
     res = 100)

PCAPlot(seurat_integrated_nuc, split.by = "sample")  + 
  ggtitle("PCA Plot by Sample_nuc seq samples only") + 
  theme(axis.text.x = element_blank()) +
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  theme(plot.title = element_text(hjust=0.5, face="bold"))

dev.off()




# Run UMAP
seurat_integrated_nuc <- RunUMAP(seurat_integrated_nuc, 
                                 dims = 1:40,
                                 reduction = "pca")

# Plot UMAP                             
DimPlot(seurat_integrated_nuc)        
Idents(seurat_integrated_nuc) <- "sample"

jpeg(filename = "figures/reduction/umap_nuc.jpeg", 
     width = 1500, height = 800,
     quality = 100,
     res = 100)

DimPlot(seurat_integrated_nuc, split.by = "sample")

dev.off()




# 3. Clustering ----
# Determine the K-nearest neighbor graph
seurat_integrated_nuc <- FindNeighbors(object = seurat_integrated_nuc, dims = 1:40)


# Determine the clusters for various resolutions                                
seurat_integrated_nuc <- FindClusters(object = seurat_integrated_nuc, 
                                      resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))


# save 
saveRDS(seurat_integrated_nuc, file = "objects/seurat_integrated_nuc_clustered.rds")
# seurat_integrated_nuc <- readRDS(file = "objects/seurat_integrated_nuc_clustered.rds")

View(names(seurat_integrated_nuc@meta.data))



# :: Clusters by resolution ----
res.list_nuc <- list(names(seurat_integrated_nuc@meta.data[17]),
                     names(seurat_integrated_nuc@meta.data[18]),
                     names(seurat_integrated_nuc@meta.data[19]),
                     names(seurat_integrated_nuc@meta.data[20]),
                     names(seurat_integrated_nuc@meta.data[21]))

res.list_nuc

# check the default assay
DefaultAssay(seurat_integrated_nuc) <- "integrated"

view_clusters_nuc <- function(res){
  Idents(seurat_integrated_nuc) <- res
  DimPlot(seurat_integrated_nuc,
          reduction = "umap",
          label = TRUE,
          label.size = 5,
          pt.size = 0.8) +
    ggtitle(paste0("Clustering of nuc seq samples at ", res)) + 
    theme(plot.title = element_text(hjust=0.5, face="bold"))
}



# re-leveling clusters
# check the number of clusters
length(unique(seurat_integrated_nuc$integrated_snn_res.0.4))    # 18
length(unique(seurat_integrated_nuc$integrated_snn_res.0.6))    # 18
length(unique(seurat_integrated_nuc$integrated_snn_res.0.8))    # 20
length(unique(seurat_integrated_nuc$integrated_snn_res.1))      # 24
length(unique(seurat_integrated_nuc$integrated_snn_res.1.4))    # 25




seurat_integrated_nuc$integrated_snn_res.0.4 <- factor(seurat_integrated_nuc$integrated_snn_res.0.4, 
                                                       levels = 0:length(unique(seurat_integrated_nuc$integrated_snn_res.0.4))-1)

seurat_integrated_nuc$integrated_snn_res.0.6 <- factor(seurat_integrated_nuc$integrated_snn_res.0.6, 
                                                       levels = 0:length(unique(seurat_integrated_nuc$integrated_snn_res.0.6))-1)


seurat_integrated_nuc$integrated_snn_res.0.8 <- factor(seurat_integrated_nuc$integrated_snn_res.0.8, 
                                                       levels = 0:length(unique(seurat_integrated_nuc$integrated_snn_res.0.8))-1)

seurat_integrated_nuc$integrated_snn_res.1 <- factor(seurat_integrated_nuc$integrated_snn_res.1, 
                                                     levels = 0:length(unique(seurat_integrated_nuc$integrated_snn_res.1))-1)

seurat_integrated_nuc$integrated_snn_res.1.4 <- factor(seurat_integrated_nuc$integrated_snn_res.1.4, 
                                                       levels = 0:length(unique(seurat_integrated_nuc$integrated_snn_res.1.4))-1)


umap.list_nuc <- list()
for(i in 1:length(res.list_nuc)){
  umap.list_nuc[[i]] <- view_clusters_nuc(res.list_nuc[[i]])
}


# export
for(i in 1:length(umap.list_nuc)){
  jpeg(filename = paste0("figures/Clustering/Clusters_umap_nuc_seq_", res.list_nuc[[i]], ".jpeg"),
       width = 1200, height = 1000,
       quality = 100,
       res = 100)
  
  print(umap.list_nuc[[i]])
  
  dev.off()
}


# :: Clustering QC ----
Idents(seurat_integrated_nuc) <- "sample"
View(table(Idents(seurat_integrated_nuc)))



library(dplyr)
n_cells_nuc.list <- list()

n_cells_nuc.list[[1]] <- FetchData(seurat_integrated_nuc, 
                                   vars = c("ident", "integrated_snn_res.0.4")) %>%
  dplyr::count(ident, integrated_snn_res.0.4) %>%
  tidyr::spread(ident, n)

n_cells_nuc.list[[2]] <- FetchData(seurat_integrated_nuc, 
                                   vars = c("ident", "integrated_snn_res.0.6")) %>%
  dplyr::count(ident, integrated_snn_res.0.6) %>%
  tidyr::spread(ident, n)


n_cells_nuc.list[[3]] <- FetchData(seurat_integrated_nuc, 
                                   vars = c("ident", "integrated_snn_res.0.8")) %>%
  dplyr::count(ident, integrated_snn_res.0.8) %>%
  tidyr::spread(ident, n)

n_cells_nuc.list[[4]] <- FetchData(seurat_integrated_nuc, 
                                   vars = c("ident", "integrated_snn_res.1")) %>%
  dplyr::count(ident, integrated_snn_res.1) %>%
  tidyr::spread(ident, n)

n_cells_nuc.list[[5]] <- FetchData(seurat_integrated_nuc, 
                                   vars = c("ident", "integrated_snn_res.1.4")) %>%
  dplyr::count(ident, integrated_snn_res.1.4) %>%
  tidyr::spread(ident, n)




# View tables
names(n_cells_nuc.list) <- res.list
View(n_cells_nuc.list[[1]])


View(table(seurat_integrated_nuc$sample))

# export
library(xlsx)
for(i in 1:length(n_cells_nuc.list)){
  write.xlsx(n_cells_nuc.list[[i]], 
             file = "results/Clustering QC/clustering_QC_ncells_by_sample_nuc_seq.xlsx", 
             sheetName = names(n_cells_nuc.list[i]),
             append = T,
             row.names = F)
}



# Explore whether clusters segregate by cell cycle phase
Idents(seurat_integrated_nuc) <- "orig.ident"
DimPlot(seurat_integrated_nuc,
        label = F, 
        split.by = "Phase",
        pt.size = 0.3)  + NoLegend()





# Segregation of clusters by various sources of uninteresting variation
# Determine metrics to plot present in seurat_integrated_nuc@meta.data
metrics <-  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "percent.mt")
DefaultAssay(seurat_integrated_nuc) <- "RNA"

cluster.qc <- function(obj, resolution){
  Idents(obj) <- resolution
  FeaturePlot(obj, 
              reduction = "umap", 
              features = metrics,
              pt.size = 0.4, 
              order = TRUE,
              min.cutoff = 'q10',             # play around with these number (min.cutoff, max.cutoff)
              # max.cutoff = 'q99',
              label = TRUE,
              label.size = 5,
              ncol = 3)
}


# export
for(i in 1:length(res.list_nuc)){
  jpeg(filename = paste0("figures/Clustering/Clustering QC metrics_nuc_seq.", res.list_nuc[[i]], ".jpeg"),
       width = 1500, height = 1000)
  
  print(cluster.qc(seurat_integrated_nuc, res.list_nuc[[i]]))
  
  dev.off()
}

