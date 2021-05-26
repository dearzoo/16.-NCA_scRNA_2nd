
# Data ready ----

split_seurat_6w <- list()
# to set different cutoffs for filtration, make a function.
filtering <- function(obj = merged_seurat, ngene, mt){
  x <- subset(x = merged_seurat,
              subset= (nFeature_RNA >= ngene) & (percent.mt <= mt))
  x <- SplitObject(x, split.by = "sample")
  x <- x[7:10]
}


split_seurat_6w[[1]] <- filtering(merged_seurat, 200, 40)
split_seurat_6w[[2]] <- filtering(merged_seurat, 250, 40)
split_seurat_6w[[3]] <- filtering(merged_seurat, 300, 40)
split_seurat_6w[[4]] <- filtering(merged_seurat, 200, 50)
split_seurat_6w[[5]] <- filtering(merged_seurat, 250, 50)
split_seurat_6w[[6]] <- filtering(merged_seurat, 300, 50)

names(split_seurat_6w) <- c("split_seurat_6w_200_40", 
                            "split_seurat_6w_250_40",
                            "split_seurat_6w_300_40",
                            "split_seurat_6w_200_50", 
                            "split_seurat_6w_250_50",
                            "split_seurat_6w_300_50")

cell_counts_fil <- list()
for(i in 1:length(split_seurat_6w)){
  cell_counts_fil[[i]] <- unlist(lapply(split_seurat_6w[[i]], ncol))
}

names(cell_counts_fil) <- names(split_seurat_6w)
View(cell_counts_fil)


get.cell.counts <- function(i){
  unlist(lapply(split_seurat_6w[[i]], ncol))
}

View(unlist(lapply(seurat_integrated_6w, ncol)))


library(purrr)
filtered_cell_counts <- map_dfr(c(1:length(split_seurat_6w)), get.cell.counts)
View(filtered_cell_counts)
rownames(filtered_cell_counts) <- names(split_seurat_6w)




# 1. Pre-process ----
# normalize and identify variable features for each dataset independently
for(i in 1:length(split_seurat_6w)){
  split_seurat_6w[[i]] <- lapply(X = split_seurat_6w[[i]], FUN = function(x) {
    x <- NormalizeData(x, verbose = TRUE)
    x <- CellCycleScoring(x, g2m.features=g2m_genes, s.features=s_genes)
    x <- SCTransform(x, vars.to.regress = c("percent.mt"))
  })
}


# 2. Integration ----
# Select the most variable features to use for integration (default: nfeatures = 2000)
integ_features_6w <- list()
for(i in 1:length(split_seurat_6w)){
  integ_features_6w[[i]] <- SelectIntegrationFeatures(object.list = split_seurat_6w[[i]], nfeatures = 3000)
}


# Prepare the SCT list object for integration
for(i in 1:length(split_seurat_6w)){
  split_seurat_6w[[i]] <- PrepSCTIntegration(object.list = split_seurat_6w[[i]],
                                       anchor.features = integ_features_6w[[i]])
}

# Find best buddies - can take a while to run
integ_anchors_6w <- list()
for(i in 1:length(split_seurat_6w)){
  integ_anchors_6w[[i]] <- FindIntegrationAnchors(object.list = split_seurat_6w[[i]],
                                                  normalization.method = "SCT",
                                                  anchor.features = integ_features_6w[[i]],
                                                  k.filter = 50,
                                                  k.score = 15)
}

# Integrate across conditions
seurat_integrated_6w <- list()
for(i in 1:length(split_seurat_6w)){
  seurat_integrated_6w[[i]] <- IntegrateData(anchorset = integ_anchors_6w[[i]],
                                             normalization.method = "SCT")
}

names(seurat_integrated_6w) <- names(split_seurat_6w)
View(seurat_integrated_6w)

# save
saveRDS(seurat_integrated_6w, file = "objects/6w only/seurat_integrated_6w.rds")
# seurat_integrated_6w <- readRDS(file = "objects/6w only/seurat_integrated_6w.rds")


# Re-leveling sample variable
names(table(seurat_integrated_6w[[1]]$sample))   # "w6_45f.02" "w6_47f.02" "w6_g10"    "w6_wt"
for(i in 1:length(seurat_integrated_6w)){
  seurat_integrated_6w[[i]]$sample <- factor(seurat_integrated_6w[[i]]$sample, 
                                             levels = c("w6_wt",
                                                        "w6_g10",
                                                        "w6_47f.02",
                                                        "w6_45f.02"))
}



# Perform integrated analysis independently
# 3. Reductions & Clustering ----
# normalize and identify variable features for each dataset independently
seurat_integrated_6w <- lapply(X = seurat_integrated_6w, FUN = function(x) {
  x <- RunPCA(x, npcs = 50, verbose = FALSE)
  x <- RunUMAP(x, reduction = "pca", dims = 1:40)
  x <- RunTSNE(x, reduction = "pca", dims = 1:40)
  x <- FindNeighbors(x, reduction = "pca", dims = 1:40)
  x <- FindClusters(x, resolution = c(0.5, 1.0))
})


# save
saveRDS(seurat_integrated_6w, file = "objects/6w only/seurat_integrated_clustered_6w.rds")
# seurat_integrated_6w <- readRDS(file = "objects/6w only/seurat_integrated_clustered_6w.rds")

names(seurat_integrated_6w[[1]]@meta.data)

# Plotting 
# Default assay setting
lapply(X=seurat_integrated_6w, FUN = function(x){
  DefaultAssay(x) <- "integrated"
})


umap_plot_6w <- lapply(X=seurat_integrated_6w, FUN = function(x){
  DimPlot(x,
          reduction = "umap",
          label = TRUE,
          label.size = 5,
          pt.size = 0.8)
})

names(umap_plot_6w) <- names(split_seurat_6w)
View(umap_plot_6w)

# exporting
for(i in 1:length(umap_plot_6w)){
  jpeg(filename = paste0("branch_6w/figures/UMAP at 1.0 of ", names(umap_plot_6w[i]), ".jpeg"), 
       res = 200, quality = 100, width = 1000, height = 1000)
  print(umap_plot_6w[[i]] + 
          ggtitle(names(umap_plot_6w[i])) + 
          theme(plot.title = element_text(hjust=0.5, face="bold"))) 
  dev.off()
}

# 4. Cell Type ----
m.mature.astro <- c("AQP4", "AGT", "GJA1")
m.fetal.astro <- c("TUBA1B")
m.oligo <- c("FGF1", "PLP1")
m.neurons <- c("DCX", "KCNJ6", "FOXP2")

marker.list <- list(m.mature.astro,
                    m.fetal.astro,
                    m.oligo,
                    m.neurons)

names(marker.list) <- c("m.mature.astro",
                        "m.fetal.astro",
                        "m.oligo",
                        "m.neurons")



# :: FeaturePlots ----
fp <- function(x, redu, markers){FeaturePlot(x,
                                             reduction = redu,
                                             features = markers, cols = c("grey", "red"),
                                             label = TRUE)
}

for (i in 1:length(seurat_integrated_6w)){
  DefaultAssay(seurat_integrated_6w[[i]]) <- "RNA"
}


# FeaturePlot_UMAP
umap.feature.list <- list()
for(i in 1:length(marker.list)){
  umap.feature.list[[i]] <- fp(seurat_integrated_6w[[1]], "umap", marker.list[[i]])
}

for(i in 1:length(marker.list)){
  umap.feature.list[[i]] <- lapply(X = seurat_integrated_6w, FUN = function(x, redu, markers){
    x <- fp(x, "umap", marker.list[[i]])
  })
}

names(umap.feature.list) <- names(marker.list)

# export as jpeg
names(umap.feature.list[[1]])
names(umap.feature.list[1])

# seurat 300_40
for(i in 1:length(umap.feature.list)){
  jpeg(filename = paste0("branch_6w/figures/cell type/seurat_300_40/FP_umap_300_40_", names(umap.feature.list[i]), ".jpeg"), 
       res = 150, quality = 100, width = 1000, height = 1000)
    print(umap.feature.list[[i]][[3]])
    dev.off()
}

jpeg(filename = paste0("branch_6w/figures/cell type/seurat_300_40/FP_umap_300_40_", names(umap.feature.list[3]), ".jpeg"), 
     res = 150, quality = 100, width = 1000, height = 600)
print(umap.feature.list[[3]][[3]])
dev.off()




tsne.feature.list <- list()
for(i in 1:length(marker.list)){
  tsne.feature.list[[i]] <- lapply(X = seurat_integrated_6w, FUN = function(x, redu, markers){
    x <- fp(x, "tsne", marker.list[[i]])
  })
}

names(tsne.feature.list) <- names(marker.list)


# export as jpeg
names(tsne.feature.list[[1]])
names(tsne.feature.list[1])

# seurat 300_40
for(i in 1:length(tsne.feature.list)){
  jpeg(filename = paste0("branch_6w/figures/cell type/seurat_300_40/FP_tsne_300_40_", names(tsne.feature.list[i]), ".jpeg"), 
       res = 150, quality = 100, width = 1000, height = 1000)
  print(tsne.feature.list[[i]][[3]])
  dev.off()
}

jpeg(filename = paste0("branch_6w/figures/cell type/seurat_300_40/FP_tsne_300_40_", names(tsne.feature.list[3]), ".jpeg"), 
     res = 150, quality = 100, width = 1000, height = 600)
print(tsne.feature.list[[3]][[3]])
dev.off()





# Monocle3 ----

# 1. CDS creation ----
# function to extract raw counts
fetch_count <- function(x){
  as.sparse(GetAssayData(x, slot = "counts"))
}

# memory
options(future.globals.maxSize = 4000 * 1024^8)

# need a count data (extract from seurat object) from split_seurat
# make a list
count.list_6w <- list()

# fetch count data from each sample
for(i in 1:length(split_seurat_6w)){
  count.list_6w[[i]] <- lapply(split_seurat_6w[[i]], fetch_count)
}

names(count.list_6w) <- paste0("count_", names(split_seurat_6w))


# check cell counts of each
for(i in 1:length(count.list_6w)){
  print(unlist(lapply(count.list_6w[[i]], ncol)))
}

# w6_wt    w6_g10 w6_47f.02 w6_45f.02 
# 114       116       181       175 
# w6_wt    w6_g10 w6_47f.02 w6_45f.02 
# 89        93       150       139 
# w6_wt    w6_g10 w6_47f.02 w6_45f.02 
# 71        81       126       120 
# w6_wt    w6_g10 w6_47f.02 w6_45f.02 
# 127       126       189       185 
# w6_wt    w6_g10 w6_47f.02 w6_45f.02 
# 94        98       155       147 
# w6_wt    w6_g10 w6_47f.02 w6_45f.02 
# 76        84       129       127 


# Cell metadata
# make a list
meta.list <- list()

for(i in 1:length(split_seurat_6w)){
  meta.list[[i]] <- lapply(X = split_seurat_6w[[i]], FUN = function(x){
    x <- x@meta.data
  })
}

names(meta.list) <- paste0("meta_", names(split_seurat_6w))
View(meta.list)


# comeback here (3/31) ====
# Gene annotation
gene.anno.list <- list()
for(i in 1:length(split_seurat_6w)){
  gene.anno.list[[i]] <- data.frame("gene_short_name" = rownames(split_seurat_6w[[i]]))
  rownames(gene.anno.list[[i]]) <- gene.anno.list[[i]][, 1]
}

names(gene.anno.list) <- paste0("gene_", names(split_seurat_6w))

View(gene.anno.list)


# Assembly 
cds_list <- list()

for(i in 1:length(count.list_6w)){
  cds_list[[i]] <- new_cell_data_set(expression_data = count.list_6w[[i]],
                                     cell_metadata = meta.list[[i]],
                                     gene_metadata = gene.anno.list[[i]])
}

names(cds_list) <- names(split_seurat_6w)
View(cds_list)


# _ Combine CDS objects ----
cds <- combine_cds(cds_list)

# 2. Pre-process ----
cds <- preprocess_cds(cds, num_dim = 100)
