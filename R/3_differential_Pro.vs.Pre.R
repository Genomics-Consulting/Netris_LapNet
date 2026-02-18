
# main --------------------------------------------------------------------

# microBulk RNAseq: human FFPE + LCM, PDAC tissues
# effect of therapy on transcriptome


# Settings ----

rm(list=ls())
gc()

suppressPackageStartupMessages({
  library(stringr)
  library(edgeR)
  library(org.Hs.eg.db)
  library(NMF)
  library(ggfortify)
  library(eulerr)
  library(pals)
  library(gridExtra)
  library(EnhancedVolcano)
  library(scales)
  library(dplyr)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(clusterProfiler)
  library(enrichplot)
  library(DOSE)
  library(pathview)
  library(ReactomePA)
  library(pathfindR)
  library(hrbrthemes)
  library(ggpubr)
  library(msigdbr)
  library(patchwork)
  
})

source("R/helper.functions.R", echo=TRUE)

set.seed(1)

out.dir <- "results/"
data.dir <- "data/"


# Post vs Pre ----

## unpaired ----

load(file="data/LapNet/DGElist_all.RData")

pdata <- y$samples
head(pdata)

# remove sample with few reads (24N0 2454-06-006)
samples.to.remove <- pdata[pdata$Structures == "Stroma", ]
samples.to.remove <- c(rownames(samples.to.remove), "06-006_Pre_Tumor")
y <- y[, -which(colnames(y) %in% samples.to.remove)]
dim(y) #  27773    28
pdata <- y$samples
head(pdata)
table(pdata$timepoint, pdata$Structures)

# Model
pdata$timepoint <- as.factor(pdata$timepoint)
pdata$timepoint <- relevel(pdata$timepoint, ref = "Pre")
design <- model.matrix(~timepoint, data = pdata)
ncol(design)
head(design)

# Filtering
keep <- filterByExpr(y, design)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y) # 13889    28

# Normalization for composition bias
y <- calcNormFactors(y)

# Estimate dispersion
y <- estimateDisp(y, design, robust = T)
head(y$samples)
plotBCV(y)

save(y, file=paste0(data.dir, "DGElist_norm.RData"))

# Differential expression
load(file=paste0(data.dir, "DGElist_norm.RData"))

pdata <- y$samples
head(pdata)
pdata$timepoint <- as.factor(pdata$timepoint)
pdata$timepoint <- relevel(pdata$timepoint, ref = "Pre")
design <- model.matrix(~timepoint, data = pdata)
ncol(design)
head(design)

fit <- glmQLFit(y, design, robust=TRUE)
qlf <- glmQLFTest(fit)
topTags(qlf)
plotQLDisp(fit)
summary(fit$df.prior)

hist(qlf$table[,"PValue"], breaks=50)
topTags(qlf)
is.de <- decideTests(qlf)
summary(is.de)
# timepointPost
# Down              28
# NotSig         13629
# Up               232

plotMD(qlf, main = "Time Point: Post vs Pre", ylim=c(-5,5))
abline(h=c(-1, 1), col="blue")

top <- topTags(qlf, sort.by = "PValue", p.value = 1, n= 50000)
top <- top$table
head(top, 20)

sort(top[top$FDR < 0.05 & top$logFC > 0, "Symbol"])
sort(top[top$FDR < 0.05 & top$logFC < 0, "Symbol"])

write.csv(top, file = paste0(out.dir, "DEGs_Post.vs.Pre.csv"))



### GSEA ----

kd <- read.csv(file = paste0(out.dir, "DEGs_Post.vs.Pre.csv"), row.names = 1)
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
head(results.table)

write.csv(results.table, file = paste0(out.dir, "GSEA_Hallmarks_Post.vs.Pre.csv"), row.names = F)
saveRDS(em2, paste0(out.dir, "GSEA_Hallmarks_Post.vs.Pre.rds"))


## paired ----

load(file="data/LapNet/DGElist_all.RData")
# keep only samples with paired Pre/Post tumor tissues
y <- y[, y$samples$Structures=="Tumor", ]
y <- y[, y$samples$ID_paired%in%c("01-007","02-001","02-003","02-014")]
dim(y) # 27773    8

pdata <- y$samples
head(pdata)
table(pdata$ID_paired, pdata$timepoint)

# Model
pdata$timepoint <- as.factor(pdata$timepoint)
pdata$timepoint <- relevel(pdata$timepoint, ref = "Pre")
design <- model.matrix(~ID_paired + timepoint, data = pdata)
ncol(design)
head(design)

# Filtering
keep <- filterByExpr(y, design)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y) # 13380    8

# Normalization for composition bias
y <- calcNormFactors(y)

# Estimate dispersion
y <- estimateDisp(y, design, robust = T)
head(y$samples)
plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
qlf <- glmQLFTest(fit, coef = "timepointPost")
topTags(qlf)
plotQLDisp(fit)
summary(fit$df.prior)

hist(qlf$table[,"PValue"], breaks=50)
topTags(qlf)
is.de <- decideTests(qlf)
summary(is.de)
plotMD(qlf, main = "Time Point: Post vs Pre", ylim=c(-5,5))
abline(h=c(-1, 1), col="blue")

top <- topTags(qlf, sort.by = "PValue", p.value = 1, n= 50000)
top <- top$table
head(top, 20)

write.csv(top, file = paste0(out.dir, "DEGs_Post.vs.Pre_paired.csv"))


### GSEA ----

kd <- read.csv(file = paste0(out.dir, "DEGs_Post.vs.Pre_paired.csv"), row.names = 1)
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

write.csv(results.table, file = paste0(out.dir, "GSEA_Hallmarks_Post.vs.Pre_paired.csv"), row.names = F)
saveRDS(em2, paste0(out.dir, "GSEA_Hallmarks_Post.vs.Pre_paired.rds"))



# end ---------------------------------------------------------------------
sessionInfo()







