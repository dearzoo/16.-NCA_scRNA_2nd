# NCA 2nd Round 6W Data Only

## Data prep

- start with merged seurat that has not been yet filtered. split it by sample and take 6w data only
- try different minimum gene counts and mt content cutoffs
  - split_seurat_6w_200_40 means minimum gene counts per cell is 200 and mito content 40% or less

**Cell Counts by Sample**
      
|                        | 6w_WT | 6w_G10 | 6w_I47F | 6w_SI |
|:----------------------:|:-----:|:------:|:-------:|:-----:|
| split_seurat_6w_200_40 |  114  |   116  |   181   |  175  |
| split_seurat_6w_250_40 |   89  |   93   |   150   |  139  |
| split_seurat_6w_300_40 |   71  |   81   |   126   |  120  |
| split_seurat_6w_200_50 |  127  |   126  |   189   |  185  |
| split_seurat_6w_250_50 |   94  |   98   |   155   |  147  |
| split_seurat_6w_300_50 |   76  |   84   |   129   |  127  |
  

## Pre-process and integration
- Normalization
- CellCycleScoring
- SCTransform
- Finding the most variable genes
- PrepSCTIntegration
- Find Integration Anchors
  - errors with default parameters
  - k value parameters were modified
- Integration across conditions (samples)

## Integrated analysis
- PCA
- UMAP
- TSNE
- Clustering (at resolution 0.5 and 1.0)
  - Not enough clusters were found
  
  * **resolution 0.5**
  
  <img src="../branch_6w/figures/UMAP at 0.5 of split_seurat_6w_200_40.jpeg" width="500">
  <img src="../branch_6w/figures/UMAP at 0.5 of split_seurat_6w_250_40.jpeg" width="500">
  <img src="../branch_6w/figures/UMAP at 0.5 of split_seurat_6w_300_40.jpeg" width="500">
  <img src="../branch_6w/figures/UMAP at 0.5 of split_seurat_6w_200_50.jpeg" width="500">
  <img src="../branch_6w/figures/UMAP at 0.5 of split_seurat_6w_250_50.jpeg" width="500">
  <img src="../branch_6w/figures/UMAP at 0.5 of split_seurat_6w_300_50.jpeg" width="500">
  
  
  * **resolution 1.0**
  
  <img src="../branch_6w/figures/UMAP at 1.0 of split_seurat_6w_200_40.jpeg" width="500">
  <img src="../branch_6w/figures/UMAP at 1.0 of split_seurat_6w_250_40.jpeg" width="500">
  <img src="../branch_6w/figures/UMAP at 1.0 of split_seurat_6w_300_40.jpeg" width="500">
  <img src="../branch_6w/figures/UMAP at 1.0 of split_seurat_6w_200_50.jpeg" width="500">
  <img src="../branch_6w/figures/UMAP at 1.0 of split_seurat_6w_250_50.jpeg" width="500">
  <img src="../branch_6w/figures/UMAP at 1.0 of split_seurat_6w_300_50.jpeg" width="500">
  


  


