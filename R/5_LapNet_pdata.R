# main --------------------------------------------------------------------

# prepare full LapNet pdata for Figures


# Settings ----

rm(list=ls())
gc()

suppressPackageStartupMessages({
  library(stringr)
  library(edgeR)
  library(dplyr)
  library(msigdbr)
  library(singscore)
})

source("R/helper.functions.R", echo=TRUE)

set.seed(1)

out.dir <- "results/"
data.dir <- "data/"



# final table ----

# this is the complete version of the normalized dataset
# with all metadata information
# only tumor samples after removing stroma and bad quality sample
load(file="data/LapNet/DGElist_norm.RData")
dim(y) # 13889    28
pdata <- y$samples
head(pdata)

logcpm <- cpm(y, log=TRUE, normalized.lib.sizes = T)
rownames(logcpm) <- y$genes$Symbol
logcpm[1:5,1:5]

## add Netrin genes ----
netrin_genes <- c("NTN1","NTN3","NTN4","NTN5",
                  "UNC5A","UNC5B","UNC5C","UNC5D",
                  "DCC","NEO1","RGMA")

setdiff(netrin_genes, rownames(logcpm)) # "NTN3" "NTN5" are missing
netrin_genes <- intersect(netrin_genes, rownames(logcpm))
netrin_genes <- sort(netrin_genes)
logcpm.sel <- logcpm[netrin_genes, ]
logcpm.sel <- as.data.frame(t(logcpm.sel))
pdata <- cbind(pdata, logcpm.sel)
head(pdata)


## add EMT Core Score ----

rankData <- rankGenes(logcpm)

# get GSEA results
gsea.results <- read.csv(file = paste0(out.dir, "GSEA_Hallmarks_Post.vs.Pre_v2.csv"), row.names = 1)
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


## add EMT Full Score ----

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


## save full object ----

pdata$timepoint <- as.factor(pdata$timepoint)
pdata$timepoint <- relevel(pdata$timepoint, ref = "Pre")

y$samples <- pdata

save(y, file=paste0(data.dir, "DGElist_norm.RData"))
write.csv(pdata, file = paste0(out.dir, "tumor.all_full.metadata.csv"), row.names = T)


# end ---------------------------------------------------------------------
sessionInfo()
