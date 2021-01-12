# Set up ----
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


# Data overview ----
# https://github.com/dearzoo/16.-NCA_scRNA_2nd


# Data read-in ----
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


# converting function
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
                      "w9_wt")              

# check dimension
unlist(lapply(data.list, nrow))
unlist(lapply(data.list, ncol))
# >>>>> Data All Read-in ====

# . ----
# Gene annotation (BioMart) ----
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
seurat.list <- lapply(good.data.list, seurat.obj, n=200)    
seurat.list <- lapply(good.data.list, seurat.obj, n=100)    
View(seurat.list)
