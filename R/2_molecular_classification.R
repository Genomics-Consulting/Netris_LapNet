
# main --------------------------------------------------------------------

# microBulk RNAseq: human FFPE + LCM, PDAC tissues
# paired fastq files

# PDACMOC was installed in a mamba environment following the instructions at:
# https://github.com/pavillos/PDACMOC
# Reference: https://www.biorxiv.org/content/10.1101/2025.03.06.641837v1.full
# check issue: https://github.com/pavillos/PDACMOC/issues/2

# PuriST function was downloaded from GitHub:
# https://github.com/lootpiz/PurIST/blob/master/R/PurIST.R
# References: https://pubmed.ncbi.nlm.nih.gov/31754050/
# https://www.jmdjournal.org/article/S1525-1578%2824%2900180-6/fulltext


# Settings ----

rm(list=ls())
gc()

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(Rsubread)
  library(Rsamtools)
  library(stringr)
  library(edgeR)
  library(org.Hs.eg.db)
  library(ggfortify)
  library(eulerr)
  library(pals)
  library(gridExtra)
  library(pals)
  library(scales)
  library(dplyr)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(hrbrthemes)
  library(ggpubr)
})

set.seed(1)

source("R/helper.functions.R", echo=TRUE)

out.dir <- "results/"
data.dir <- "data/"


# Molecular classifiers for PDAC ----

## prepare TCGA PAAD ----

# downloading PAAD to get count level
query <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
PAAD_data <- GDCprepare(query)
counts_matrix <- assay(PAAD_data)
print(dim(counts_matrix))
print(counts_matrix[1:5, 1:5])
sample_metadata <- colData(PAAD_data)
gene_info <- rowData(PAAD_data)
print(head(gene_info))
print(head(sample_metadata))
colnames(sample_metadata)
table(sample_metadata$`paper_mRNA Moffitt clusters (76 High Purity Samples Only)  1basal  2classical`)
table(sample_metadata$`paper_mRNA Moffitt clusters (All 150 Samples) 1basal  2classical`)
table(sample_metadata$`paper_mRNA Collisson clusters (All 150 Samples) 1classical 2exocrine 3QM`)
table(sample_metadata$`paper_mRNA Bailey Clusters (All 150 Samples) 1squamous 2immunogenic 3progenitor 4ADEX`)

# saving for later use
save(counts_matrix, sample_metadata, gene_info, 
     file= paste0(data.dir, "TCGA_PAAD_counts_metadata.RData"))

all(rownames(counts_matrix) == gene_info$gene_id) # TRUE
rownames(counts_matrix) <- gene_info$gene_name
save(counts_matrix, sample_metadata, gene_info, 
     file= paste0(data.dir, "TCGA_PAAD_counts.RData"))


## PDACMOC ------------------------------------------------------------

****************************************
  shell script in mamba environment
****************************************
  
mamba activate PDACMOC

R

library(PDACMOC)

### test on example data ----

file1 <- system.file('examples', 'example.R', package = 'PDACMOC')
dir <- file.path(dirname(file1), '../saved_workspaces/example.RData')
#load(dir)
rm(file1)
# Look for the PDACMOC environment automatically
reticulate::use_condaenv("PDACMOC", required = TRUE)
file2 <- system.file('training_data', 'all_datasets_corrected.csv', package = 'PDACMOC')
# in the example 'samples' is in counts format with Ensembl IDs as row names
samples <- read.csv(file2, row.names = 1, check.names = FALSE)
rm(file2)

new_samples <- PDACMOC:::import.and.normalize(samples, batch = FALSE, gene_id = 'EnsemblID')
# this changed EnsemblIDs to gene symbols
classification_tumor <- PDACMOC:::omni.classify(samples, batch = FALSE, gene_id = 'EnsemblID',
                                                classifier = c('Collisson', 'Moffitt', 'Bailey',
                                                               'Puleo', 'Chan-Seng-Yue', 'PDAConsensus'))

classification_all <- PDACMOC:::omni.classify(samples, batch = FALSE, gene_id = 'EnsemblID',
                                              classifier = c('Collisson', 'Moffitt', 'Bailey',
                                                             'Puleo', 'Chan-Seng-Yue', 'PDAConsensus'),
                                              stroma = TRUE,
                                              stroma_classifier = c('Moffitt',
                                                                    'Maurer',
                                                                    'PDAConsensus'))
q()

mamba deactivate


### test on TCGA-PAAD ----

mamba activate PDACMOC
R
library(PDACMOC)
setwd('')
load("TCGA_PAAD_counts.RData")
reticulate::use_condaenv("PDACMOC", required = TRUE)

classification_tumor <- PDACMOC:::omni.classify(counts_matrix, batch = FALSE, gene_id = 'GeneSymbol',
                                                classifier = c('Collisson', 'Moffitt', 'Bailey',
                                                               'Puleo', 'Chan-Seng-Yue', 'PDAConsensus'))
save(classification_tumor, file="TCGA_PAAD_classification_tumor.RData")

q()
mamba deactivate

# check here

load(paste0(data.dir, "TCGA_PAAD_counts_metadata.RData"))
head(sample_metadata)
table(sample_metadata$`paper_mRNA Moffitt clusters (76 High Purity Samples Only)  1basal  2classical`)
table(sample_metadata$`paper_mRNA Moffitt clusters (All 150 Samples) 1basal  2classical`)
table(sample_metadata$`paper_mRNA Collisson clusters (All 150 Samples) 1classical 2exocrine 3QM`)
table(sample_metadata$`paper_mRNA Bailey Clusters (All 150 Samples) 1squamous 2immunogenic 3progenitor 4ADEX`)

load(paste0(data.dir, "TCGA_PAAD_classification_tumor.RData"))
names(classification_tumor$Tumor)
#"Collisson"     "Moffitt"       "Bailey"        "Puleo"         "Chan-Seng-Yue" "PDAConsensus"

head(classification_tumor$Tumor$Collisson)
lapply(classification_tumor$Tumor, dim)

# compare results

table(classification_tumor$Tumor$Collisson$Predicted.subtype, 
      sample_metadata$`paper_mRNA Collisson clusters (All 150 Samples) 1classical 2exocrine 3QM`)
#               1  2  3
# Classical     48  6 28
# QM             0  0  5
# Exocrine-like  6 56  1

table(classification_tumor$Tumor$Moffitt$Predicted.subtype, 
      sample_metadata$`paper_mRNA Moffitt clusters (All 150 Samples) 1basal  2classical`)
#            1  2
# Classical   5 72
# Basal-like 60 13

table(classification_tumor$Tumor$Bailey$Predicted.subtype, 
      sample_metadata$`paper_mRNA Bailey Clusters (All 150 Samples) 1squamous 2immunogenic 3progenitor 4ADEX`)
#              1  2  3  4
# Squamous     1  0  0  0
# Progenitor  15  4 29  1
# Immunogenic 11 20 11  2
# ADEX         4  4 13 35

# this is similar to what was found by the authors (see issue above)


### test in Linehan-seq ----

mamba activate PDACMOC
R
library(PDACMOC)

load("DGElist_Linehan.all.RData")

counts <- y$counts
rownames(counts) <- y$genes$Symbol
counts[1:5,1:5]
dim(counts)

reticulate::use_condaenv("PDACMOC", required = TRUE)

classification_tumor <- PDACMOC:::omni.classify(counts, batch = FALSE, gene_id = 'GeneSymbol',
                                                classifier = c('Collisson', 'Moffitt', 'Bailey',
                                                               'Puleo', 'Chan-Seng-Yue', 'PDAConsensus'))
save(classification_tumor, file="Linehan_classification_tumor.RData")

q()
mamba deactivate


# check here (from script #4)

load("DGElist_Linehan.all.RData")
pdata <- y$samples
head(pdata)

load(paste0("Linehan_classification_tumor.RData"))
names(classification_tumor$Tumor)
#"Collisson"     "Moffitt"       "Bailey"        "Puleo"         "Chan-Seng-Yue" "PDAConsensus"
head(classification_tumor$Tumor$Collisson)
head(classification_tumor$Tumor$PDAConsensus)

lapply(classification_tumor$Tumor, dim)

# compare results

table(classification_tumor$Tumor$Collisson$Predicted.subtype, 
      pdata$Collisson)
#                 Classical Exocrine-like QM-PDA
#Classical            18             2     20
#QM                    0             0      0
#Exocrine-like         7            18      1

table(classification_tumor$Tumor$Moffitt$Predicted.subtype, 
      pdata$Moffitt)
#               Basal-like Classical
#Classical           0        56
#Basal-like          6         4

table(classification_tumor$Tumor$Bailey$Predicted.subtype, 
      pdata$Bailey)
#             ADEX Immunogenic Pancreatic Progenitor Squamous
# Squamous       0           0                     0        1
# Progenitor     3          24                    14       14
# Immunogenic    0           0                     0        0
# ADEX           7           1                     1        1

# this is overall similar to what was found by the authors
# although no Immunogenic samples were detected with Bailey



### run in LapNet ----

mamba activate PDACMOC
R
library(PDACMOC)

load("DGElist_all.RData")

counts <- y$counts
rownames(counts) <- y$genes$Symbol
counts[1:5,1:5]
dim(counts) # 27773    35

reticulate::use_condaenv("PDACMOC", required = TRUE)

classification_tumor <- PDACMOC:::omni.classify(counts, batch = FALSE, gene_id = 'GeneSymbol',
                                                classifier = c('Collisson', 'Moffitt', 'Bailey',
                                                               'Puleo', 'Chan-Seng-Yue', 'PDAConsensus'),
                                                stroma = TRUE,
                                                stroma_classifier = c('Moffitt',
                                                                      'Maurer',
                                                                      'PDAConsensus')
                                                )
save(classification_tumor, file="LapNet_classification_tumor.RData")

q()
mamba deactivate


## Purist ----

# initially used a version of PurIST downloaded from GitHub
# that uses binary penalization. It resulted in a single basal-like sample only.
# next switched to a version adapted to work with continuous data (log2 ratio)
# that requires normalized (non-logarightmic) data
# this resulted in 3 basal-like samples (including the one obtained with the binary version)
# tested in TCGA with ~30 samples classified as basal-like (below)

### validate on TCGA-PAAD ----

load(paste0(data.dir, "TCGA_PAAD_counts.RData"))
dim(counts_matrix)
head(counts_matrix[,1:5])
rownames(counts_matrix) <- gene_info$gene_name

res.tcga <- apply(counts_matrix, 2, PurIST_LogRatio)
#FALSE  TRUE 
#162    21
res.tcga <- apply(counts_matrix, 2, PurIST_binary)
table(res.tcga > 0.5)
#FALSE  TRUE 
#156    27

load(paste0(data.dir, "TCGA_PAAD_counts_metadata.RData"))
head(sample_metadata)

all(names(res.tcga) == rownames(sample_metadata)) # TRUE

table(res.tcga > 0.5, sample_metadata$`paper_mRNA Moffitt clusters (76 High Purity Samples Only)  1basal  2classical`)
table(res.tcga > 0.5, sample_metadata$`paper_mRNA Moffitt clusters (All 150 Samples) 1basal  2classical`)
table(res.tcga > 0.5, sample_metadata$`paper_mRNA Collisson clusters (All 150 Samples) 1classical 2exocrine 3QM`)
table(res.tcga > 0.5, sample_metadata$`paper_mRNA Bailey Clusters (All 150 Samples) 1squamous 2immunogenic 3progenitor 4ADEX`)

# only partially matches Moffitt classification

### validate on MetaGxPancreas ----

library(MetaGxPancreas)

pancreasData <- loadPancreasDatasets()
names(pancreasData)

SEs <- pancreasData$SummarizedExperiments
names(SEs)
lapply(SEs, dim)
# how many total samples?
sum(unlist(lapply(SEs, ncol))) # 1710 samples

# check for purist genes in each dataset
purist.biomarkers = c("GPR87", "REG4", 
               "KRT6A", "ANXA10",
               "BCAR3", "GATA6",
               "PTGES", "CLDN18",
               "ITGA3", "LGALS4",
               "C16orf74", "DDC", 
               "S100A2", "SLC40A1", 
               "KRT5", "CLRN3")
lapply(SEs, function(SE) {
  exprs <- assays(SE)$exprs
  setdiff(purist.biomarkers, rownames(exprs))
})

# BALAGURANATH_SumExp and HAMIDI_SumExp are missing several biomarkers
#names(SEs)
#SEs <- SEs[-c(2, 16)]

# apply purist to each dataset where all genes are present

res <- PurIST_Batch_Predict(assays(SEs[[1]])$exprs)
table(res$PurIST_Probability > 0.5)
hist(res$PurIST_Probability, breaks = 20)
table(res$Subtype_Call)

purist.results <- lapply(SEs, function(SE) {
  exprs <- assays(SE)$exprs
  res <- PurIST_Batch_Predict(exprs)
})

# check results
names(purist.results)
lapply(purist.results, function(res) {
  table(res$Subtype_Call)
})
lapply(purist.results, function(res) {
  table(res$PurIST_Probability > 0.5)
})


# which proportion of samples are basal-like in each dataset?

purist.proportions <- sapply(purist.results, function(res) {
  if (is.null(res)) {
    return(NA)
  } else {
    return(mean(res$PurIST_Probability > 0.5))
  }
})
summary(purist.proportions)
# 2% of samples are basal like


# update pdata ---------------------------------------------------------------

load(file=paste0("data/LapNet/DGElist_all.RData"))
pdata <- y$samples
head(pdata)
counts <- y$counts
rownames(counts) <- y$genes$Symbol
counts[1:5,1:5]
dim(counts)

load(paste0(data.dir, "LapNet_classification_tumor.RData"))

head(classification_tumor$Proportions)

names(classification_tumor$Tumor)
#"Collisson"     "Moffitt"       "Bailey"        "Puleo"         "Chan-Seng-Yue" "PDAConsensus"
head(classification_tumor$Tumor$Collisson)
lapply(classification_tumor$Tumor, dim)
names(classification_tumor$Stroma)
# "Moffitt"      "Maurer"       "PDAConsensus"
head(classification_tumor$Stroma$PDAConsensus)
lapply(classification_tumor$Stroma, dim)

pdata$Collisson <- classification_tumor$Tumor$Collisson$Predicted.subtype
pdata$Moffitt <- classification_tumor$Tumor$Moffitt$Predicted.subtype
pdata$Bailey <- classification_tumor$Tumor$Bailey$Predicted.subtype
pdata$Puleo <- classification_tumor$Tumor$Puleo$Predicted.subtype
pdata$ChanSengYue <- classification_tumor$Tumor$`Chan-Seng-Yue`$Predicted.subtype
pdata$PDAConsensus <- classification_tumor$Tumor$PDAConsensus$Predicted.subtype
pdata$Stroma.Moffitt <- classification_tumor$Stroma$Moffitt$Predicted.subtype
pdata$Stroma.Maurer <- classification_tumor$Stroma$Maurer$Predicted.subtype
pdata$Stroma.Consensus <- classification_tumor$Stroma$PDAConsensus$Predicted.subtype

pdata$Collisson_probability <- classification_tumor$Tumor$Collisson$Probability
pdata$Moffitt_probability <- classification_tumor$Tumor$Moffitt$Probability
pdata$Bailey_probability <- classification_tumor$Tumor$Bailey$Probability
pdata$Puleo_probability <- classification_tumor$Tumor$Puleo$Probability
pdata$ChanSengYue_probability <- classification_tumor$Tumor$`Chan-Seng-Yue`$Probability
pdata$PDAConsensus_probability <- classification_tumor$Tumor$PDAConsensus$Probability
pdata$Stroma.Moffitt_probability <- classification_tumor$Stroma$Moffitt$Probability
pdata$Stroma.Maurer_probability <- classification_tumor$Stroma$Maurer$Probability
pdata$Stroma.Consensus_probability <- classification_tumor$Stroma$PDAConsensus$Probability

pdata <- cbind(pdata, classification_tumor$Proportions)
head(pdata)


## save updated pdata
y$samples <- pdata
save(y, file=paste0(data.dir, "/DGElist_all.RData"))

write.csv(pdata, file = paste0(out.dir, "LapNet.all_molecular_classification.csv"), row.names = T)


## ESTIMATE tumor purity estimation ----

library(tidyestimate)

load(file=paste0(data.dir, "DGElist_all.RData"))
pdata <- y$samples
head(pdata)

counts <- cpm(y, log=F, normalized.lib.sizes = T)
rownames(counts) <- y$genes$Symbol
counts[1:5,1:5]
dim(counts)

filtered <- filter_common_genes(counts, 
                                id = "hgnc_symbol", 
                                tidy = FALSE, 
                                tell_missing = TRUE, 
                                find_alias = TRUE)
scored <- estimate_score(filtered,
                         is_affymetrix = TRUE)
head(scored)
plot_purity(scored, is_affymetrix = TRUE)

pdata <- cbind(pdata, scored)

y$samples <- pdata
save(y, file=paste0(data.dir, "DGElist_all.RData"))

write.csv(pdata, file = paste0(out.dir, "LapNet.all_molecular_classification.csv"), row.names = T)


# similar results with rpkm
rpkm_counts <- rpkm(y, log=F, normalized.lib.sizes = T, gene.length = y$genes$Length)
rownames(rpkm_counts) <- y$genes$Symbol
rpkm_counts[1:5,1:5]
dim(rpkm_counts)

filtered <- filter_common_genes(rpkm_counts, 
                                id = "hgnc_symbol", 
                                tidy = FALSE, 
                                tell_missing = TRUE, 
                                find_alias = TRUE)
scored <- estimate_score(filtered,
                         is_affymetrix = TRUE)
head(scored)
plot_purity(scored, is_affymetrix = TRUE)


## PurIST ----

load(file=paste0(data.dir, "DGElist_all.RData"))
pdata <- y$samples
head(pdata)

counts <- cpm(y, log=F, normalized.lib.sizes = T)
rownames(counts) <- y$genes$Symbol
counts[1:5,1:5]
dim(counts)
# similar results with rpkm
#rpkm_counts <- rpkm(y, log=F, normalized.lib.sizes = T, gene.length = y$genes$Length)
#rownames(rpkm_counts) <- y$genes$Symbol
#rpkm_counts[1:5,1:5]
#dim(rpkm_counts)

res <- apply(counts, 2, PurIST_LogRatio)
head(res)
#res <- apply(rpkm_counts, 2, PurIST_LogRatio)
# if the PBP is greater than 0.5, the tumor subtype is determined to be a basal- like subtype and 
# if the PBP if less than or equal to 0.5, the tumor subtype is determined to be a classical subtype.
table(res > 0.5)
hist(res, breaks = 20)

pdata$PurIST.score_LogRatio <- res
pdata$PurIST_LogRatio <- ifelse(res > 0.5, "Basal-like", "Classical")
head(pdata)

table(pdata$PurIST_LogRatio)
# 3 basal-like samples now (but one of them stromal)
y$samples <- pdata
save(y, file=paste0(data.dir, "DGElist_all.RData"))



# end ---------------------------------------------------------------------
sessionInfo()




