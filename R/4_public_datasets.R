# main --------------------------------------------------------------------

# GSE253260
# Nicolle et al: https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(24)00409-2/fulltext
# Predictive genomic and transcriptomic analysis on endoscopic ultrasound-guided fine needle aspiration materials 
# from primary pancreatic adenocarcinoma: a prospective multicentre study
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253260
# processed from raw counts

# GSE131050
# Linehan dataset described in the PurIST paper (https://pubmed.ncbi.nlm.nih.gov/31754050/)
# corresponds to Timothy M Nywening et al: https://pubmed.ncbi.nlm.nih.gov/27055731/
# Phase 1b study targeting tumour associated macrophages with CCR2 inhibition plus 
# FOLFIRINOX in locally advanced and borderline resectable pancreatic cancer
# fine needle aspiration (FNA) tumor biopsies were collected at baseline and after completion of treatment cycle 2
#	https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131050
#	processed from fastq files (BioProject	PRJNA542388) using Rsubread and featureCounts


# Settings ----

rm(list=ls())
gc()

suppressPackageStartupMessages({
  library(GEOquery)
  library(stringr)
  library(edgeR)
  library(org.Hs.eg.db)
  library(pals)
  library(gridExtra)
  library(scales)
  library(dplyr)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(hrbrthemes)
  library(ggpubr)
  library(survival)
  library(survminer)
  library(singscore)
  library(msigdbr)
})

source("R/helper.functions.R", echo=TRUE)

set.seed(1)

out.dir <- "results/"
data.dir <- "data/"


# GSE253260 ----

## Prepare data ----

gse <- getGEO('GSE253260',GSEMatrix=TRUE)

show(gse) # 397 samples
pdata <- pData(gse$GSE253260_series_matrix.txt.gz)
colnames(pdata)
head(pdata)
pdata <- pdata[, c(1,2,53:64)]
colnames(pdata) <- c("title", "geo_accession", "stage", "ecog", "firsttreatment",
                     "historicpatient_id", "normpatient_id", "OS_time", "OS_event",
                     "PFS_time", "PFS_event", "sample_id", "Sex", "sizeprimary") 

table(pdata$stage)
# borderline locAdvanced        meta  resectable 
# 55         170         105          61 
table(pdata$firsttreatment)
table(pdata$OS_event)
pdata$OS_time <- as.numeric(as.character(pdata$OS_time))
pdata$OS_event <- as.numeric(as.character(pdata$OS_event))
pdata$PFS_time <- as.numeric(as.character(pdata$PFS_time))
pdata$PFS_event <- as.numeric(as.character(pdata$PFS_event))

# downloaded from GEO
raw.counts <- read.delim(paste0(data.dir, "GSE253260_BACAP.rawct.tsv.gz"), row.names = 1)
dim(raw.counts) # 60671 396
head(raw.counts)[,1:5]

setdiff(colnames(raw.counts), pdata$title)
setdiff(pdata$title, colnames(raw.counts)) # BPC137T

pdata <- pdata[pdata$title %in% colnames(raw.counts), ]

all(colnames(raw.counts) == pdata$title) # TRUE

# get gene symbols
raw.counts$Symbol <- mapIds(org.Hs.eg.db, rownames(raw.counts), keytype="ENSEMBL", column="SYMBOL")
raw.counts <- raw.counts[!is.na(raw.counts$Symbol), ]
dim(raw.counts) # 36688 397
raw.counts <- raw.counts[!duplicated(raw.counts$Symbol), ]
dim(raw.counts) # 36594 397
rownames(raw.counts) <- raw.counts$Symbol
raw.counts$Symbol <- NULL

boxplot(log2(raw.counts[,1:50]+1), las = 2)

y <- DGEList(counts= raw.counts, samples = pdata)

save(y, file=paste0(data.dir, "DGElist_GSE253260.RData"))


## Process subset ----

# select LAPC samples treated with Folfirinox

head(pdata)
rownames(pdata) <- pdata$title
pdata <- pdata[pdata$stage %in% c("locAdvanced"), ]
dim(pdata) # 170 14
raw.counts <- raw.counts[, rownames(pdata)]
all(colnames(raw.counts) == pdata$title) # TRUE

y <- DGEList(counts= raw.counts, samples = pdata)

save(y, file=paste0(data.dir, "DGElist_GSE253260_LAPC.RData"))

# only Folfirinox
table(pdata$firsttreatment)
table(is.na(pdata$firsttreatment))
samples.to.keep <- pdata[!is.na(pdata$firsttreatment), ]
samples.to.keep <- samples.to.keep[samples.to.keep$firsttreatment == "FOLFIRINOX", ] # n=43
samples.to.keep <- rownames(samples.to.keep)
y <- y[, which(colnames(y) %in% samples.to.keep)]

dim(y) # 36594   43
pdata <- y$samples
table(pdata$stage)
table(pdata$stage, pdata$firsttreatment)

# survival data is given in months
pdata$OS_time <- as.numeric(as.character(pdata$OS_time)) * 30
pdata$PFS_time <- as.numeric(as.character(pdata$PFS_time)) * 30

# normalization
y <- calcNormFactors(y)


## EMT ----

logcpm <- cpm(y, log=TRUE, normalized.lib.sizes = T)
all(rownames(pdata)==colnames(logcpm)) # TRUE

rankData <- rankGenes(logcpm)
dim(rankData)
rankData[1:15,1:5]

## Core EMT

gsea.results <- read.csv(file = paste0(out.dir, "GSEA_Hallmarks_Post.vs.Pre.csv"), row.names = 1)
head(gsea.results)
core.genes <- gsea.results["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "core_enrichment"]
core.genes <- unlist(strsplit(core.genes, "/"))
sort(core.genes) # 45 genes
setdiff(core.genes, rownames(logcpm)) # all genes present
intersect(core.genes, rownames(logcpm))

scoredf <- simpleScore(rankData, upSet = core.genes)
head(scoredf)
pdata$core.EMT_score <- scoredf$TotalScore
head(pdata)

# full EMT 

library(escape)
species = "Homo sapiens"
gene.set.list <- list(hallmark = c("H"))
gene.sets <- getGeneSets(species = species, 
                         library = gene.set.list[[1]][1], subcategory = gene.set.list[[1]][2])
sel.gene.sets <- gene.sets[c("HALLMARK-EPITHELIAL-MESENCHYMAL-TRANSITION")]
names(sel.gene.sets)
signature.genes <- as.vector(sel.gene.sets$`HALLMARK-EPITHELIAL-MESENCHYMAL-TRANSITION`) # 200
scoredf <- simpleScore(rankData, upSet = signature.genes)
head(scoredf)
pdata$singscore_EMT <- scoredf$TotalScore


## Netrin genes ----

netrin_genes <- c("NTN1","NTN3","NTN4","NTN5",
                  "UNC5A","UNC5B","UNC5C","UNC5D",
                  "DCC","NEO1","RGMA","ADORA2B")
intersect(netrin_genes, rownames(logcpm)) # 12/12 present 

netrin_exp <- logcpm[rownames(logcpm) %in% netrin_genes, ]
netrin_exp <- as.data.frame(t(netrin_exp))
colnames(netrin_exp) <- paste0("exp_", colnames(netrin_exp))
head(netrin_exp)
all(rownames(netrin_exp) == pdata$title)
pdata <- cbind(pdata, netrin_exp)
head(pdata)


y$samples <- pdata
save(y, file="data/LapNet/DGElist_GSE253260_LAPC_full.pdata.RData")



# GSE131050 -----------------------------------------------------

## Prepare Data ----

gds <- getGEO("GSE131050")
show(gds)
show(pData(phenoData(gds[[1]]))) # Linehan_seq (n=66)
pdata0 <- pData(phenoData(gds[[1]]))[, 1:2]
head(pdata0)

# pdata from the PurIST paper supplementary data
pdata1 <- read.csv(paste0(data.dir, "pdata_Linehan.csv"))
head(pdata1)
all(pdata1$ID==pdata0$title)
pdata1$accession <- rownames(pdata0)

# metadata from SRA website
pdata2 <- read.csv(paste0(data.dir, "Linehan_SraRunTable.csv"))
head(pdata2)
rownames(pdata2) <- pdata2$GEO_Accession..exp.
pdata2 <- pdata2[pdata1$accession, ]

all(rownames(pdata2)==pdata1$accession)
pdata1$Run <- pdata2$Run

write.csv(pdata1, file=paste0(data.dir, "pdata_Linehan_full.csv"), row.names = F)


## counts after Rsubread + featureCounts
load("data/LapNet/counts_Linehan.RData")
head(fc$counts)
counts <- fc$counts
colnames(counts) <- gsub(".bam","",colnames(counts))

pdata <- read.csv(paste0(data.dir, "pdata_Linehan_full.csv"))
head(pdata)
rownames(pdata) <- pdata$Run

all(colnames(counts) == rownames(pdata)) # TRUE

table(pdata$Treatment)
# FOLFIRINOX FOLFIRINOX+PF04136309 
# 14                    52 
table(pdata$Treatment, pdata$Pre.Post)
#                         Post Pre
# FOLFIRINOX               5   9
# FOLFIRINOX+PF04136309   24  28


## Process edgeR ----

y <- DGEList(counts= counts, samples = pdata)
dim(y) # 28395    66

# Add gene annotation
y$genes <- fc$annotation[, c("Length","GeneID"), drop=FALSE]
y$genes$Symbol <- mapIds(org.Hs.eg.db, rownames(y), keytype="ENTREZID", column="SYMBOL")
y$genes$EntrezID <- rownames(y)
head(y$genes)
y <- y[!is.na(y$genes$Symbol), ]
dim(y) # 27773    66

save(y, file= paste0(data.dir, "DGElist_Linehan.all.RData"))


## Process paired ----

# select only Folfirinox paired samples
paired.samples <- c("S1124.02","S1124.03","S1124.07","S1124.08","S1124.13")
sel.samples <- rownames(pdata[pdata$Patient_ID %in% paired.samples, ]) # 10 samples

y <- y[, sel.samples]
pdata <- y$samples
head(pdata)
table(pdata$Pre.Post, pdata$Patient_ID) # 5 matched patients

# Model
pdata$Pre.Post <- as.factor(pdata$Pre.Post)
pdata$Pre.Post <- relevel(pdata$Pre.Post, ref = "Pre")
design <- model.matrix(~Patient_ID + Pre.Post, data = pdata)
ncol(design)
head(design)

# Filtering
keep <- filterByExpr(y, design)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y) # 17592

# Normalization for composition bias
y <- calcNormFactors(y)
colors <- rainbow(4)
points <- 15:18
group <- as.numeric(as.factor(pdata$Pre.Post))
plotMDS.DGEList(y, pch = points[group], col= colors[group], main= "MDS Plot: normalized data", cex = 2)
legend("topleft", legend=levels(group), pch=points, col=colors, ncol=1)

# Estimate dispersion
y <- estimateDisp(y,design, robust = T)
y$samples
plotBCV(y)

save(y, file= paste0(data.dir, "DGElist_Linehan.Folfirinox_paired.RData"))


## Netrin genes ----

load(paste0(data.dir, "DGElist_Linehan.Folfirinox_paired.RData"))
dim(y) # 17592    10
pdata <- y$samples
head(pdata)

logcpm <- cpm(y, log=TRUE, normalized.lib.sizes = T)
rownames(logcpm) <- y$genes$Symbol
logcpm[1:5,1:5]

netrin_genes <- c("NTN1","NTN3","NTN4","NTN5",
                  "UNC5A","UNC5B","UNC5C","UNC5D",
                  "DCC","NEO1","RGMA","ADORA2B")

setdiff(netrin_genes, rownames(logcpm)) # "DCC"
netrin_genes <- intersect(netrin_genes, rownames(logcpm))
netrin_genes <- sort(netrin_genes)
logcpm.sel <- logcpm[netrin_genes, ]
logcpm.sel <- as.data.frame(t(logcpm.sel))
pdata <- cbind(pdata, logcpm.sel)
head(pdata)


## EMT ----
library(singscore)
rankData <- rankGenes(logcpm)

# get GSEA results from LapNet Pre vs Post analysis
gsea.results <- read.csv(file = paste0("results/LapNet/GSEA_Hallmarks_Post.vs.Pre_v2.csv"), row.names = 1)
head(gsea.results)
core.genes <- gsea.results["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "core_enrichment"]
core.genes <- unlist(strsplit(core.genes, "/"))
sort(core.genes) # 45 genes

setdiff(core.genes, rownames(logcpm))
intersect(core.genes, rownames(logcpm))

# singscore
scoredf <- simpleScore(rankData, upSet = core.genes)
head(scoredf)
pdata$core.EMT_score <- scoredf$TotalScore

head(pdata)


## add EMT Full Score

rankData <- rankGenes(logcpm)

library(escape)

# hallmarks
species = "Homo sapiens"
gene.set.list <- list(hallmark = c("H"))
gene.sets <- getGeneSets(species = species, 
                         library = gene.set.list[[1]][1], subcategory = gene.set.list[[1]][2])
names(gene.sets)
sel.gene.sets <- gene.sets[c("HALLMARK-EPITHELIAL-MESENCHYMAL-TRANSITION")]
names(sel.gene.sets)
signature.genes <- as.vector(sel.gene.sets$`HALLMARK-EPITHELIAL-MESENCHYMAL-TRANSITION`) # 200
scoredf <- simpleScore(rankData, upSet = signature.genes)
head(scoredf)
pdata$singscore_EMT <- scoredf$TotalScore
head(pdata)

y$samples <- pdata

save(y, file= paste0(data.dir, "DGElist_Linehan.Folfirinox_paired.RData"))
write.csv(pdata, file = paste0(out.dir, "Linehan_folfirinox_paired_metadata.csv"), row.names = T)




## Integration -------------------------------------------------------------

load(file=paste0("data/LapNet/DGElist_Linehan.all.RData"))
Linehan <- y
dim(Linehan)
load(file=paste0(data.dir, "DGElist_all.RData")) # LapNet
Lapnet <- y
dim(Lapnet)
rm(y)

pdata.Linehan <- Linehan$samples
head(pdata.Linehan)
# select only Folfirinox samples
sel.samples <- rownames(pdata.Linehan[pdata.Linehan$Treatment=="FOLFIRINOX", ]) # 14
Linehan <- Linehan[, sel.samples]
pdata.Linehan <- Linehan$samples
head(pdata.Linehan)
table(pdata.Linehan$Pre.Post)

pdata.Lapnet <- Lapnet$samples
head(pdata.Lapnet)
# remove stromal samples and sample with few reads (24N0 2454-06-006)
samples.to.remove <- pdata.Lapnet[pdata.Lapnet$Structures == "Stroma", ]
samples.to.remove <- c(rownames(samples.to.remove), "06-006_Pre_Tumor")
Lapnet <- Lapnet[, -which(colnames(Lapnet) %in% samples.to.remove)]
pdata.Lapnet <- Lapnet$samples
head(pdata.Lapnet)
table(pdata.Lapnet$timepoint)

counts.Linehan <- Linehan$counts
dim(counts.Linehan) # 27773   14
counts.Linehan[1:5,1:10]
counts.Lapnet <- Lapnet$counts
dim(counts.Lapnet) # 27773   28
counts.Lapnet[1:5,1:10]

# keep only common genes
common.genes <- intersect(rownames(counts.Linehan), rownames(counts.Lapnet)) # 20000
counts.Linehan <- counts.Linehan[common.genes, ]
counts.Lapnet <- counts.Lapnet[common.genes, ]
counts <- cbind(counts.Linehan, counts.Lapnet)
head(counts)

# combined pdata
pdata.Linehan$study <- "Linehan_seq"
pdata.Linehan <- pdata.Linehan[, c("Treatment","Pre.Post","Patient_ID","study")]
colnames(pdata.Linehan) <- c("Treatment","Timepoint","Patient_ID","Study")
head(pdata.Linehan)

pdata.Lapnet$study <- "LapNet_seq"
pdata.Lapnet$Treatment <- "FOLFIRINOX_NP137"
pdata.Lapnet <- pdata.Lapnet[, c("Treatment","timepoint","ID_paired","study")]
colnames(pdata.Lapnet) <- c("Treatment","Timepoint","Patient_ID", "Study")
head(pdata.Lapnet)

pdata <- rbind(pdata.Linehan, pdata.Lapnet)
head(pdata)
table(pdata$Study)
table(pdata$Timepoint, pdata$Patient_ID)
table(pdata$Treatment, pdata$Timepoint)
# Post Pre
# FOLFIRINOX          5   9
# FOLFIRINOX_NP137    6  22

all(rownames(pdata) == colnames(counts)) # TRUE

y <- DGEList(counts= counts, samples = pdata, group = as.factor(pdata$Treatment))

# Add gene annotation
y$genes <- as.data.frame(rownames(y$counts))
y$genes$Symbol <- mapIds(org.Hs.eg.db, rownames(y), keytype="ENTREZID", column="SYMBOL")
colnames(y$genes) <- c("EntrezID","Symbol")
head(y$genes)
y <- y[!is.na(y$genes$Symbol), ]
dim(y) # 27773 42

save(y, file=paste0(data.dir, "DGElist_Linehan.Lapnet_integrated.RData"))


## Folfirinox.NP137 vs Folfirinox ----

# Model
group <- y$samples$group
batch <- y$samples$Study
design <- model.matrix(~0+group+batch)
head(design)

# Filtering
keep <- filterByExpr(y, design)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y) # 18696

# Normalization for composition bias
y <- calcNormFactors(y)

# Estimate dispersion
y <- estimateDisp(y,design, robust = T)
plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)
plotQLDisp(fit)
summary(fit$df.prior)

contrast <- makeContrasts(groupFOLFIRINOX_NP137-groupFOLFIRINOX, levels=design)
res <- glmQLFTest(fit, contrast= contrast)
hist(res$table[,"PValue"], breaks=50)
topTags(res)
is.de <- decideTests(res)
summary(is.de)
par(mfrow=c(1,2), mar=c(8,6,4,4))
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"), 
       main = "FOLFIRINOX_NP137 vs FOLFIRINOX\nMD plot",
       legend="topright", pch = 20)

top <- topTags(res, sort.by = "PValue", p.value = 1, n= 20000)
# data for GSEA:
write.csv(top, file = paste0(out.dir, "integrated_FOLFIRINOX_NP137.vs.FOLFIRINOX_unpaired.csv"))


### GSEA ----

kd <- read.csv(paste0(out.dir, "integrated_FOLFIRINOX_NP137.vs.FOLFIRINOX_unpaired.csv"), row.names = 1)
head(kd)
# considering p value and fold change for ranking the genes
kd$neg.logP <- -log10(kd$PValue)
kd$gsea <- kd$neg.logP * kd$logFC
kd <- kd[order(kd$gsea, decreasing = T), ]
rownames(kd) <- kd$EntrezID
geneList = kd[,"gsea"]
names(geneList) = as.character(rownames(kd))
head(geneList)
tail(geneList)

# Get MSigDB Hallmark gene sets
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

em2 <- GSEA(geneList, 
            pvalueCutoff = 1,
            seed = 4,
            TERM2GENE = m_t2g)
head(em2, 10)
edox <- setReadable(em2, 'org.Hs.eg.db', 'ENTREZID')
results.table <- as.data.frame(edox)

write.csv(results.table, file = paste0(out.dir, "GSEA_Hallmarks_Folfi.NP137.vs.Folfi.csv"), row.names = F)
saveRDS(em2, paste0(out.dir, "GSEA_Hallmarks_Folfi.NP137.vs.Folfi.rds"))


## Treatment-specific effect ----

load(file=paste0(data.dir, "DGElist_Linehan.Lapnet_integrated.RData"))
dim(y)

# select only paired samples
pdata <- y$samples
head(pdata)
table(pdata$Patient_ID, pdata$Timepoint)
duplicated.patients <- names(which(table(pdata$Patient_ID)==2))
sel.samples <- rownames(pdata[pdata$Patient_ID %in% duplicated.patients, ])
y <- y[, sel.samples]
dim(y)

# The model strategy is the Blocked Design Approach (~ 0 + group2 with Patient_ID as a block). 
# as this is a two-treatment, two-timepoint paired (nested) study 
# where each patient is only in one treatment arm (data collinearity issues).

pdata <- y$samples
pdata$Patient_ID <- factor(pdata$Patient_ID)
pdata$Timepoint <- factor(pdata$Timepoint)
pdata$Treatment <- factor(pdata$Treatment)
table(pdata$Treatment, pdata$Timepoint, pdata$Patient_ID)
pdata$group2 <- factor(paste(pdata$Timepoint, pdata$Treatment, sep = "."))
head(pdata)
table(pdata$group2)

design <- model.matrix(~ 0 + group2, pdata)

head(design)
colnames(design)
keep <- filterByExpr(y, design)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y) # 17155    18
y <- calcNormFactors(y)
y <- estimateDisp(y,design, robust = T, block = Patient_ID)
plotBCV(y)

save(y, file=paste0(data.dir, "DGElist_Linehan.Lapnet_integrated_paired.RData"))

fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)
plotQLDisp(fit)
summary(fit$df.prior)

# Define the contrast for the difference in the timepoint effect (Interaction)
# This compares the magnitude of the Post vs Pre change between the two treatments, 
# isolating the unique effect of NP137 that goes beyond the effect of FOLFIRINOX.
contrast_Interaction <- makeContrasts((group2Post.FOLFIRINOX_NP137 - group2Pre.FOLFIRINOX_NP137) - 
                                        (group2Post.FOLFIRINOX - group2Pre.FOLFIRINOX), 
                                      levels = design)
qlf_Interaction <- glmQLFTest(fit, contrast = contrast_Interaction)
results_Interaction <- topTags(qlf_Interaction, n = Inf)$table
head(results_Interaction, 50)
hist(results_Interaction[,"PValue"], breaks=50)
is.de <- decideTests(qlf_Interaction)
summary(is.de)
par(mfrow=c(1,2), mar=c(8,6,4,4))
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"), 
       main = "FOLFIRINOX_NP137 vs FOLFIRINOX\nMD plot",
       legend="topright", pch = 20)

# data for GSEA:
write.csv(results_Interaction, file = paste0(out.dir, "integrated_FOLFIRINOX_NP137.vs.FOLFIRINOX_paired.interaction.csv"))


### GSEA ----

kd <- read.csv(paste0(out.dir, "integrated_FOLFIRINOX_NP137.vs.FOLFIRINOX_paired.interaction.csv"), row.names = 1)
head(kd)
# considering p value and fold change for ranking the genes
kd$neg.logP <- -log10(kd$PValue)
kd$gsea <- kd$neg.logP * kd$logFC
kd <- kd[order(kd$gsea, decreasing = T), ]
rownames(kd) <- kd$EntrezID
geneList = kd[,"gsea"]
names(geneList) = as.character(rownames(kd))
head(geneList)
tail(geneList)

# Get MSigDB Hallmark gene sets
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

em2 <- GSEA(geneList, 
            pvalueCutoff = 1,
            seed = 4,
            TERM2GENE = m_t2g)
head(em2, 20)
edox <- setReadable(em2, 'org.Hs.eg.db', 'ENTREZID')
results.table <- as.data.frame(edox)

write.csv(results.table, file = paste0(out.dir, "GSEA_Hallmarks_Folfi.NP137.vs.Folfi_paired.interaction.csv"), row.names = F)
saveRDS(em2, paste0(out.dir, "GSEA_Hallmarks_Folfi.NP137.vs.Folfi_paired.interaction.rds"))



# end ---------------------------------------------------------------------
sessionInfo()

