# 1. Set up ----
library(Seurat)
library(cowplot)
library(umap)
library(ggplot2)
library(reticulate)
library(Matrix)
library(data.table)
library(dplyr)
library(SeuratData)
library(Rtsne)
library(readr)
library(stringr)
library(Matrix.utils)


# 2. Data overview ----
# https://github.com/dearzoo/16.-NCA_scRNA_2nd


# 2.1. Data read-in ----
read.in <- function(x){read.table(x, as.is = T, header = F)[, 1]}

# sample 01: data file invalid
w6_45f.01 <- readMM("data/Raw Data/NCA_3_6/01_6w_45f.01_NCA_3_6-CTCTCTAT/NCA_3_6-CTCTCTAT.mtx")
dim(w6_45f.01)

# sample 02: data file invalid
w6_47f.01 <- readMM("data/Raw Data/NCA_3_6/02_6w_47f.01_NCA_3_6-TATCCTCT/NCA_3_6-TATCCTCT.mtx")
dim(w6_47f.01)

# sample 03
w6_g10 <- readMM("data/Raw Data/NCA_3_6/03_6w_g10_NCA_3_6-GTAAGGAG/NCA_3_6-GTAAGGAG.mtx")
dim(w6_g10)    # 58333 24100
w6_g10.cn <- read.in("data/Raw Data/NCA_3_6/03_6w_g10_NCA_3_6-GTAAGGAG/NCA_3_6-GTAAGGAG.mtx.colnames")
w6_g10.rn <- read.in("data/Raw Data/NCA_3_6/03_6w_g10_NCA_3_6-GTAAGGAG/NCA_3_6-GTAAGGAG.mtx.rownames")

# sample 04
w6_wt <- readMM("data/Raw Data/NCA_3_6/04_6w_wt_NCA_3_6-ACTGCATA/NCA_3_6-ACTGCATA.mtx")
dim(w6_wt)   # 58333 16362
w6_wt.cn <- read.in("data/Raw Data/NCA_3_6/04_6w_wt_NCA_3_6-ACTGCATA/NCA_3_6-ACTGCATA.mtx.colnames")
w6_wt.rn <- read.in("data/Raw Data/NCA_3_6/04_6w_wt_NCA_3_6-ACTGCATA/NCA_3_6-ACTGCATA.mtx.rownames")

# sample 05
w6_45f.02 <- readMM("data/Raw Data/NCA_3_6/05_6w_45f.02_NCA_3_6-AAGGAGTA/NCA_3_6-AAGGAGTA.mtx")
dim(w6_45f.02)     # 58333 14932
w6_45f.02.cn <- read.in("data/Raw Data/NCA_3_6/05_6w_45f.02_NCA_3_6-AAGGAGTA/NCA_3_6-AAGGAGTA.mtx.colnames")
w6_45f.02.rn <- read.in("data/Raw Data/NCA_3_6/05_6w_45f.02_NCA_3_6-AAGGAGTA/NCA_3_6-AAGGAGTA.mtx.rownames")

# sample 06
w6_47f.02 <- readMM("data/Raw Data/NCA_3_6/06_6w_47f.02_NCA_3_6-CTAAGCCT/NCA_3_6-CTAAGCCT.mtx")
dim(w6_47f.02)    # 58333 13586
w6_47f.02.cn <- read.in("data/Raw Data/NCA_3_6/06_6w_47f.02_NCA_3_6-CTAAGCCT/NCA_3_6-CTAAGCCT.mtx.colnames")
w6_47f.02.rn <- read.in("data/Raw Data/NCA_3_6/06_6w_47f.02_NCA_3_6-CTAAGCCT/NCA_3_6-CTAAGCCT.mtx.rownames")

# sample 07
w3_45f <- readMM("data/Raw Data/NCA_3_6/07_3w_45f_NCA_3_6-CGTCTAAT/NCA_3_6-CGTCTAAT.mtx")
dim(w3_45f)    # 58333 17538
w3_45f.cn <- read.in("data/Raw Data/NCA_3_6/07_3w_45f_NCA_3_6-CGTCTAAT/NCA_3_6-CGTCTAAT.mtx.colnames")
w3_45f.rn <- read.in("data/Raw Data/NCA_3_6/07_3w_45f_NCA_3_6-CGTCTAAT/NCA_3_6-CGTCTAAT.mtx.rownames")

# sample 08
w3_47f <- readMM("data/Raw Data/NCA_3_6/08_3w_47f_NCA_3_6-TCTCTCCG/NCA_3_6-TCTCTCCG.mtx")
dim(w3_47f)     # 58333 24271
w3_47f.cn <- read.in("data/Raw Data/NCA_3_6/08_3w_47f_NCA_3_6-TCTCTCCG/NCA_3_6-TCTCTCCG.mtx.colnames")
w3_47f.rn <- read.in("data/Raw Data/NCA_3_6/08_3w_47f_NCA_3_6-TCTCTCCG/NCA_3_6-TCTCTCCG.mtx.rownames")

# sample 09
w3_g10 <- readMM("data/Raw Data/NCA_3_6/09_3w_g10_NCA_3_6-TCGACTAG/NCA_3_6-TCGACTAG.mtx")
dim(w3_g10)     # 58333 27762
w3_g10.cn <- read.in("data/Raw Data/NCA_3_6/09_3w_g10_NCA_3_6-TCGACTAG/NCA_3_6-TCGACTAG.mtx.colnames")
w3_g10.rn <- read.in("data/Raw Data/NCA_3_6/09_3w_g10_NCA_3_6-TCGACTAG/NCA_3_6-TCGACTAG.mtx.rownames")

# sample 10
w3_wt <- readMM("data/Raw Data/NCA_3_6/10_3w_wt_NCA_3_6-TTCTAGCT/NCA_3_6-TTCTAGCT.mtx")
dim(w3_wt)     # 58333 28987
w3_wt.cn <- read.in("data/Raw Data/NCA_3_6/10_3w_wt_NCA_3_6-TTCTAGCT/NCA_3_6-TTCTAGCT.mtx.colnames")
w3_wt.rn <- read.in("data/Raw Data/NCA_3_6/10_3w_wt_NCA_3_6-TTCTAGCT/NCA_3_6-TTCTAGCT.mtx.rownames")

# sample 11
w3_45f_nuc <- readMM("data/Raw Data/NCA_3_6/11_3w_45f_nuc_NCA_3_6-CCTAGAGT/NCA_3_6-CCTAGAGT.mtx")
dim(w3_45f_nuc)   # 58333 30112
w3_45f_nuc.cn <- read.in("data/Raw Data/NCA_3_6/11_3w_45f_nuc_NCA_3_6-CCTAGAGT/NCA_3_6-CCTAGAGT.mtx.colnames")
w3_45f_nuc.rn <- read.in("data/Raw Data/NCA_3_6/11_3w_45f_nuc_NCA_3_6-CCTAGAGT/NCA_3_6-CCTAGAGT.mtx.rownames")

# sample 12
w3_47f_nuc <- readMM("data/Raw Data/NCA_3_6/12_3w_47f_nuc_NCA_3_6-GCGTAAGA/NCA_3_6-GCGTAAGA.mtx")
dim(w3_47f_nuc)    # 58333 26540
w3_47f_nuc.cn <- read.in("data/Raw Data/NCA_3_6/12_3w_47f_nuc_NCA_3_6-GCGTAAGA/NCA_3_6-GCGTAAGA.mtx.colnames")
w3_47f_nuc.rn <- read.in("data/Raw Data/NCA_3_6/12_3w_47f_nuc_NCA_3_6-GCGTAAGA/NCA_3_6-GCGTAAGA.mtx.rownames")



# 9 and 12 weeks samples: need to find complementary and reverse seq
barcodes <- list("CTATTAAG",        # sample13
                 "AAGGCTAT",        # sample14
                 "GAGCCTTA",        # sample15
                 "TTATGCGA",        # sample16
                 "GGAGGTAA",        # sample17
                 "CATAACTG",        # sample18
                 "AGTAAAGG",        # sample19
                 "TCCGTCTC")        # sample21

# reverse the barcords
library(tcR)
barcodes_rv <- lapply(barcodes, reverse.string)
unlist(barcodes)
unlist(barcodes_rv)


# converting function (complementary)
convert <- function(x){
  chartr("ATGC","TACG",x)
}


# switch each letter to the mapped complementary base
barcodes_rv_complementary <- lapply(barcodes_rv, convert)
unlist(barcodes_rv_complementary)

# read in sample 13 through sample 21
# sample 13
w12_45f <- readMM("data/Raw Data/NCA_9_12Weeks/13_NCA_9_12Weeks-CTTAATAG/NCA_9_12Weeks-CTTAATAG.mtx")
dim(w12_45f)     # 58333 22041
w12_45f.cn <- read.in("data/Raw Data/NCA_9_12Weeks/13_NCA_9_12Weeks-CTTAATAG/NCA_9_12Weeks-CTTAATAG.mtx.colnames")
w12_45f.rn <- read.in("data/Raw Data/NCA_9_12Weeks/13_NCA_9_12Weeks-CTTAATAG/NCA_9_12Weeks-CTTAATAG.mtx.rownames")

# sample 14
w12_47f <- readMM("data/Raw Data/NCA_9_12Weeks/14_NCA_9_12Weeks-ATAGCCTT/NCA_9_12Weeks-ATAGCCTT.mtx")
dim(w12_47f)    # 58333 21814
w12_47f.cn <- read.in("data/Raw Data/NCA_9_12Weeks/14_NCA_9_12Weeks-ATAGCCTT/NCA_9_12Weeks-ATAGCCTT.mtx.colnames")
w12_47f.rn <- read.in("data/Raw Data/NCA_9_12Weeks/14_NCA_9_12Weeks-ATAGCCTT/NCA_9_12Weeks-ATAGCCTT.mtx.rownames")

# sample 15
w12_g10 <- readMM("data/Raw Data/NCA_9_12Weeks/15_NCA_9_12Weeks-TAAGGCTC/NCA_9_12Weeks-TAAGGCTC.mtx")
dim(w12_g10)    # 58333    26
w12_g10.cn <- read.in("data/Raw Data/NCA_9_12Weeks/15_NCA_9_12Weeks-TAAGGCTC/NCA_9_12Weeks-TAAGGCTC.mtx.colnames")
w12_g10.rn <- read.in("data/Raw Data/NCA_9_12Weeks/15_NCA_9_12Weeks-TAAGGCTC/NCA_9_12Weeks-TAAGGCTC.mtx.rownames")

# sample 16
w12_wt <- readMM("data/Raw Data/NCA_9_12Weeks/16_NCA_9_12Weeks-TCGCATAA/NCA_9_12Weeks-TCGCATAA.mtx")
dim(w12_wt)     # 58333 18201
w12_wt.cn <- read.in("data/Raw Data/NCA_9_12Weeks/16_NCA_9_12Weeks-TCGCATAA/NCA_9_12Weeks-TCGCATAA.mtx.colnames")
w12_wt.rn <- read.in("data/Raw Data/NCA_9_12Weeks/16_NCA_9_12Weeks-TCGCATAA/NCA_9_12Weeks-TCGCATAA.mtx.rownames")

# sample 17
w9_45f <- readMM("data/Raw Data/NCA_9_12Weeks/17_NCA_9_12Weeks-TTACCTCC/NCA_9_12Weeks-TTACCTCC.mtx")
dim(w9_45f)   # 58333     6
w9_45f.cn <- read.in("data/Raw Data/NCA_9_12Weeks/17_NCA_9_12Weeks-TTACCTCC/NCA_9_12Weeks-TTACCTCC.mtx.colnames") 
w9_45f.rn <- read.in("data/Raw Data/NCA_9_12Weeks/17_NCA_9_12Weeks-TTACCTCC/NCA_9_12Weeks-TTACCTCC.mtx.rownames") 

# sample 18
w9_47f <- readMM("data/Raw Data/NCA_9_12Weeks/18_NCA_9_12Weeks-CAGTTATG/NCA_9_12Weeks-CAGTTATG.mtx")
dim(w9_47f)   # 58333 19521
w9_47f.cn <- read.in("data/Raw Data/NCA_9_12Weeks/18_NCA_9_12Weeks-CAGTTATG/NCA_9_12Weeks-CAGTTATG.mtx.colnames")
w9_47f.rn <- read.in("data/Raw Data/NCA_9_12Weeks/18_NCA_9_12Weeks-CAGTTATG/NCA_9_12Weeks-CAGTTATG.mtx.rownames")

# sample 19
w9_g10 <- readMM("data/Raw Data/NCA_9_12Weeks/19_NCA_9_12Weeks-CCTTTACT/NCA_9_12Weeks-CCTTTACT.mtx")
dim(w9_g10)    # 58333 18992
w9_g10.cn <- read.in("data/Raw Data/NCA_9_12Weeks/19_NCA_9_12Weeks-CCTTTACT/NCA_9_12Weeks-CCTTTACT.mtx.colnames")
w9_g10.rn <- read.in("data/Raw Data/NCA_9_12Weeks/19_NCA_9_12Weeks-CCTTTACT/NCA_9_12Weeks-CCTTTACT.mtx.rownames")

# sample 21
w9_wt <- readMM("data/Raw Data/NCA_9_12Weeks/21_NCA_9_12Weeks-GAGACGGA/NCA_9_12Weeks-GAGACGGA.mtx")
dim(w9_wt)    # 58333 19955
w9_wt.cn <- read.in("data/Raw Data/NCA_9_12Weeks/21_NCA_9_12Weeks-GAGACGGA/NCA_9_12Weeks-GAGACGGA.mtx.colnames")
w9_wt.rn <- read.in("data/Raw Data/NCA_9_12Weeks/21_NCA_9_12Weeks-GAGACGGA/NCA_9_12Weeks-GAGACGGA.mtx.rownames")


# make a list of datasets (sample03~sample12, sample01 & 02 are invalid)
data.list <- list(w6_g10,                  # sample 03
                  w6_wt,                   # sample 04
                  w6_45f.02,               # sample 05
                  w6_47f.02,               # sample 06
                  w3_45f,                  # sample 07
                  w3_47f,                  # sample 08
                  w3_g10,                  # sample 09
                  w3_wt,                   # sample 10
                  w3_45f_nuc,              # sample 11
                  w3_47f_nuc,              # sample 12
                  w12_45f,                 # sample 13
                  w12_47f,                 # sample 14
                  w12_g10,                 # sample 15
                  w12_wt,                  # sample 16
                  w9_45f,                  # sample 17
                  w9_47f,                  # sample 18
                  w9_g10,                  # sample 19
                  w9_wt)                   # sample 21
                  
   

# name list elements
names(data.list) <- c("w6_g10",              
                      "w6_wt",               
                      "w6_45f.02",               
                      "w6_47f.02",
                      "w3_45f",
                      "w3_47f",                  
                      "w3_g10",                  
                      "w3_wt",                   
                      "w3_45f_nuc",              
                      "w3_47f_nuc",
                      "w12_45f",                 # sample 13
                      "w12_47f",                 # sample 14
                      "w12_g10",                 # sample 15
                      "w12_wt",                  # sample 16
                      "w9_45f",                  # sample 17
                      "w9_47f",                  # sample 18
                      "w9_g10",                  # sample 19
                      "w9_wt")                   # sample 21

# check dimension
unlist(lapply(data.list, nrow))
View(unlist(lapply(data.list, ncol)))
# >>>>> Data All Read-in ====

# . ----
# 2.2. Gene annotation (BioMart) ----
# before doing BioMart, check all rownames to see they are all identical

# check if all rn vectors are identical
rn.list <- list(w6_g10.rn,                  # sample 03
                w6_wt.rn,                   # sample 04
                w6_45f.02.rn,               # sample 05
                w6_47f.02.rn,               # sample 06
                w3_45f.rn,                  # sample 07
                w3_47f.rn,                  # sample 08
                w3_g10.rn,                  # sample 09
                w3_wt.rn,                   # sample 10
                w3_45f_nuc.rn,              # sample 11
                w3_47f_nuc.rn,              # sample 12
                w12_45f.rn,                 # sample 13
                w12_47f.rn,                 # sample 14
                w12_g10.rn,                 # sample 15
                w12_wt.rn,                  # sample 16
                w9_45f.rn,                  # sample 17
                w9_47f.rn,                  # sample 18
                w9_g10.rn,                  # sample 19
                w9_wt.rn)                   # sample 21


for (i in 1:(length(rn.list)-1)){
  print(identical(rn.list[[i]], rn.list[[i+1]]))
}
# result: All TRUE


# getBM
# Optimize memory space before processing and split the object by sample for iterative process
options(future.globals.maxSize = 4000 * 1024^2)

# Retrieve BioMart
ensembl <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",host = "www.ensembl.org")

# hgnc symbols matching with rownames
bm <- biomaRt::getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id",
                     values = w6_g10.rn,
                     mart= ensembl)

head(bm)

# BioMart cleaning
# 01. empty symbol
ids <- which(nchar(bm$hgnc_symbol)==0)
length(ids)     # 19382
bm <- bm[-ids,]

# 02. duplicated
ids.dup <- which(duplicated(bm$hgnc_symbol))
length(ids.dup)     # 12
bm <- bm[-ids.dup,]

# 03. final indices of ensembl ids (good.ids will be used for all rownames)
good.ids <- match(bm$ensembl_gene_id, w6_g10.rn)

# check
identical(bm$ensembl_gene_id, w6_g10.rn[good.ids])    # TRUE


# _ Assembly ----
## 01. Select only rows that are good.ids in each dataset
## 02. Assemble rownames: rownames are hgnc column of the cleaned bm
## 03. Assemble colnames


# Assembly function
# take only good.ids
cleaning <- function(x){
  x <- x[good.ids, ]
}

# make a good.data.list out of data.list
good.data.list <- list()
good.data.list <- lapply(data.list, cleaning)  # good.data.list is the list of data with clean row names (no NAs, no duplicates in hgnc symbol)

# plug in rownames (BioMart symbol) to the matrices
for(i in 1:length(good.data.list)){
  rownames(good.data.list[[i]]) <- bm$hgnc_symbol
}

# colnames
cn.list <- list(w6_g10.cn,                  # sample 03
                w6_wt.cn,                   # sample 04
                w6_45f.02.cn,               # sample 05
                w6_47f.02.cn,               # sample 06
                w3_45f.cn,                  # sample 07
                w3_47f.cn,                  # sample 08
                w3_g10.cn,                  # sample 09
                w3_wt.cn,                   # sample 10
                w3_45f_nuc.cn,              # sample 11
                w3_47f_nuc.cn,              # sample 12
                w12_45f.cn,                 # sample 13
                w12_47f.cn,                 # sample 14
                w12_g10.cn,                 # sample 15
                w12_wt.cn,                  # sample 16
                w9_45f.cn,                  # sample 17
                w9_47f.cn,                  # sample 18
                w9_g10.cn,                  # sample 19
                w9_wt.cn)                   # sample 21

# plug in col names
for(i in 1:length(good.data.list)){
  colnames(good.data.list[[i]]) <- cn.list[[i]]
}


# Check: no change in column numbers after cleaning. cleaning was only for rows.
identical(unlist(lapply(good.data.list, ncol)), unlist(lapply(data.list, ncol)))


# 3. Seurat object ----
seurat.obj <- function(x, n){CreateSeuratObject(counts = x, 
                                                project = "NCA_2R", 
                                                min.cells = 10, min.features = n)}

# make Seurat objects
library(Seurat)
seurat.list <- list()
seurat.list <- lapply(good.data.list, seurat.obj, n=200) # selected for downstream
seurat.list <- lapply(good.data.list, seurat.obj, n=100) 
seurat.list <- lapply(good.data.list, seurat.obj, n=50)  
View(seurat.list)


# cell counts 
unlist(lapply(seurat.list, ncol))
View(unlist(lapply(seurat.list, ncol)))

# working with min.features = 100
# add metadata (to be edited)
for(i in 1:length(seurat.list)){
  seurat.list[[i]]$sample <- names(seurat.list[i])
  seurat.list[[i]]$age <- "6w"
  seurat.list[[i]]$line <- "g10"
  seurat.list[[i]]$date <- "8/4/2020"
  seurat.list[[i]]$space <- "whole"
  seurat.list[[i]]$tech <- "inDrop"
}

# edit age
for(i in 5:10){
  seurat.list[[i]]$age <- "3w"
}

for(i in 11:14){
  seurat.list[[i]]$age <- "12w"
}

for(i in 15:18){
  seurat.list[[i]]$age <- "9w"
}


# edit line
for(i in c(2, 8, 14, 18)){
  seurat.list[[i]]$line <- "wt"
}


for(i in c(3, 5, 9, 11, 15)){
  seurat.list[[i]]$line <- "45f"
}

for(i in c(4, 6, 10, 12, 16)){
  seurat.list[[i]]$line <- "47f"
}


# edit date
for(i in 5:10){
  seurat.list[[i]]$date <- "8/21/2020"
}

for(i in 11:14){
  seurat.list[[i]]$date <- "10/02/2020"
}

for(i in 15:18){
  seurat.list[[i]]$date <- "10/08/2020"
}


# edit space
for (i in 9:10){
  seurat.list[[i]]$space <- "nuc-seq"
}


# Add content to metadata
for(i in 1:length(seurat.list)){
  seurat.list[[i]]$percent.mt <- PercentageFeatureSet(seurat.list[[i]], pattern = "^MT-")
}

# check if the metadata is correct
View(seurat.list)

for (i in 1:length(seurat.list)){
  print(table(seurat.list[[i]]@meta.data$age))
}

for (i in 1:length(seurat.list)){
  print(table(seurat.list[[i]]@meta.data$space))
}

for (i in 1:length(seurat.list)){
  print(table(seurat.list[[i]]@meta.data$date))
}

for (i in 1:length(seurat.list)){
  print(table(seurat.list[[i]]@meta.data$line))
}

# Re-order the seurat.list
View(names(seurat.list))
seurat.list_ordered <- seurat.list[c(8, 7, 6, 5, 10, 9,            # 3 weeks
                                     2, 1, 4, 3,                   # 6 weeks
                                     18, 17, 16, 15,               # 9 weeks
                                     14, 13, 12, 11)]              # 12 weeks

View(seurat.list_ordered)

# check
all(names(seurat.list) %in% names(seurat.list_ordered))
all(unlist(lapply(seurat.list, nrow)) %in% unlist(lapply(seurat.list_ordered, nrow)))
all(unlist(lapply(seurat.list, ncol)) %in% unlist(lapply(seurat.list_ordered, ncol)))


# all_set. replace
seurat.list <- seurat.list_ordered
rm(seurat.list_ordered)

# VlnPlot Function setting
violin <- function(seurat){
  VlnPlot(seurat,
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
          ncol = 3, 
          pt.size = 0)
}


pdf(file = "figures/violin/qc_violin.pdf")
for(i in c(1:13, 15:length(seurat.list))){
  print(violin(seurat.list[[i]]))
}
dev.off()

pdf(file = "figures/violin/qc_violin_with_names.pdf")
for(i in c(1:13, 15:length(seurat.list))){
  print(violin(seurat.list[[i]])+ggtitle(names(seurat.list[i])))
}
dev.off()


# . ----
# 4. QC ----
# Merge (w9_45 will be excluded because its seurat object dimension is 0*3)
View(names(seurat.list))                                        
merged_seurat <- merge(x = seurat.list[[1]],
                       y = list(seurat.list[[2]],
                                seurat.list[[3]],
                                seurat.list[[4]],
                                seurat.list[[5]],
                                seurat.list[[6]],
                                seurat.list[[7]],
                                seurat.list[[8]],
                                seurat.list[[9]],
                                seurat.list[[10]],
                                seurat.list[[11]],
                                seurat.list[[12]],
                                seurat.list[[13]],
                                # seurat.list[[14]],    # w9_45F
                                seurat.list[[15]],
                                seurat.list[[16]],     # w12_g10
                                seurat.list[[17]],
                                seurat.list[[18]]),
                       add.cell.id = c(names(seurat.list[1]), names(seurat.list[c(2:13, 15:18)])))





# save merged seurat
saveRDS(object = merged_seurat, file = "objects/merged_seurat.rds")
# merged_seurat <- readRDS(file = "objects/merged_seurat.rds")             
  

# check
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)

length(unique(merged_seurat$sample))
table(merged_seurat$sample)

# sample re-leveling
class(merged_seurat$sample)      # character

merged_seurat$sample <- factor(merged_seurat$sample,
                               levels = c("w3_wt", "w3_g10", "w3_47f", "w3_45f", "w3_47f_nuc", "w3_45f_nuc",
                                          "w6_wt", "w6_g10", "w6_47f.02", "w6_45f.02",
                                          "w9_wt", "w9_g10", "w9_47f",
                                          "w12_wt", "w12_g10", "w12_47f", "w12_45f"))

table(merged_seurat$sample)
length(unique(merged_seurat$sample))

# _ QC metrics set up ----
# Add number of genes per UMI
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)


# if I want to add more factors for QC metrics and work on a separate metadata,
metadata <- merged_seurat@meta.data

# add new columns to the metadata
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
head(metadata)


# _ QC visualization ----

# Cell Counts by Sample
library(ggplot2)
library(tibble)


jpeg(filename = "figures/QC/QC metrics/ncells.jpeg", 
     width = 800, height = 600,
     quality = 100,
     pointsize = 12,
     res = 100)
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
dev.off()

# Number of UMIs/transcripts per cell
jpeg(filename = "figures/QC/QC metrics/umi_counts.jpeg", 
     width = 800, height = 600,
     quality = 100,
     pointsize = 12,
     res = 100)

metadata %>% 
  ggplot(aes(color=sample, x=nCount_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 300)
dev.off()

# Detected genes per cell: distribution
jpeg(filename = "figures/QC/QC metrics/unique_genes_density.jpeg", 
     width = 800, height = 600,
     quality = 100,
     pointsize = 12,
     res = 100)
metadata %>% 
  ggplot(aes(color=sample, x=nFeature_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = c(30, 50, 100)) +
  ggtitle("Unique Genes per Cell") +
  theme(plot.title = element_text(hjust=0.5, face="bold"))
dev.off()

# Detected genes per cell: distribution via boxplot
## linear
jpeg(filename = "figures/QC/QC metrics/unique_genes_box.jpeg", 
     width = 1000, height = 700,
     quality = 100,
     pointsize = 12,
     res = 100)
metadata %>% 
  ggplot(aes(x=sample, y=nFeature_RNA, fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nFeature_RNAs per Cell") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
dev.off()

## log
jpeg(filename = "figures/QC/QC metrics/unique_genes_box_log.jpeg", 
     width = 1000, height = 700,
     quality = 100,
     pointsize = 12,
     res = 100)
metadata %>% 
  ggplot(aes(x=sample, y=log10(nFeature_RNA), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nFeature_RNAs per Cell_Log Scaled") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
dev.off()


# Number of genes per UMI
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
jpeg(filename = "figures/QC/QC metrics/umi_genes_mt.jpeg", 
     width = 1200, height = 1000,
     quality = 100,
     pointsize = 12,
     res = 100)
metadata %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 300) +
  geom_hline(yintercept = 200) +
  facet_wrap(~sample)
dev.off()

# Complexity of the gene expression: Detected genes per UMI
jpeg(filename = "figures/QC/QC metrics/complexity.jpeg", 
     width = 800, height = 600,
     quality = 100,
     pointsize = 12,
     res = 100)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill= sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) +
  ggtitle("Complexity") +
  theme(plot.title = element_text(hjust=0.5, face="bold"))
dev.off()


# unique gene number (summary statistics)
summary.stat <- list()
for(i in 1:length(seurat.list)){
  summary.stat[[i]] <- summary(seurat.list[[i]]@meta.data$nFeature_RNA)
}

names(summary.stat) <- names(seurat.list)
lapply(summary.stat, print)

med_gene_count <- list()
for(i in 1:length(seurat.list)){
  med_gene_count[[i]] <- median(seurat.list[[i]]@meta.data$nFeature_RNA)
}


# comeback here ====
# order of the table is weird
names(med_gene_count) <- names(seurat.list)
class(unlist(med_gene_count))
med.gene.count_table <- data.frame("median_gene_counts" = unlist(med_gene_count))

View(med.gene.count_table)

median_gene_count <- data.frame("sample" = names(seurat.list), "med_unique_gene" = unlist(med_gene_count))
View(median_gene_count)



# comeback ====

# :: Cell level filtration ----

# 1. Sample exclusion: w12_g10 to be excluded
length(table(merged_seurat$sample))
unlist(lapply(seurat.list, ncol))

# 2. Filter out low quality reads using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat,
                            subset= (nFeature_RNA >= 200) &
                              (nFeature_RNA <= 4000) &
                              (percent.mt <= 40))


saveRDS(filtered_seurat, file = "objects/Joseph/filtered_seurat.rds")
filtered_seurat <- readRDS(file = "objects/Joseph/filtered_seurat.rds")


# 3. QC 2nd phase ----
# _ Crude process for unwanted variations ----
## simple linear regression against a cell cycle score
# Rough normalization of the counts because raw counts are not comparable
seurat_phase <- NormalizeData(filtered_seurat)


# # Cell cycle effect
# # Load cell cycle markers (Gained at Harvard Core Training)
# load("../99. Materials/Training/Harvard scRNA-seq Analysis/single_cell_rnaseq/data/cycle.rda")

# Score cells for cell cycle (rename it because it's only to evaluate cell cycle effects)
seurat_phase <- CellCycleScoring(seurat_phase, 
                                   g2m.features = g2m_genes, 
                                   s.features = s_genes)

# View cell cycle scores and phases assigned to cells                             
View(seurat_phase@meta.data)  
dim(seurat_phase@meta.data)  # [1] 2377   12

# After scoring the cells for cell cycle, we would like to determine whether cell cycle is a major source of variation in our dataset using PCA. 
# 1. choose the most variable features
# 2. scale the data

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                       selection.method = "vst",
                                       nfeatures = 2000, 
                                       verbose = FALSE)


