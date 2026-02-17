# main --------------------------------------------------------------------

# Figure 3b: LapNet Post vs Pre, volcano
# Figure 3c: LapNet Post vs Pre, GSEA dotplot + ridgeplot
# Figure 3d: LapNet Post vs Pre, EMT violins
# Figure 3e: Linehan Post vs Pre, EMT violins
# Figure 3f: Integrated LapNet/Linehan, GSEA dotplot + ridgeplot
# Figure 3g: LapNet correlation EMT:PFS

# Figure 4a: Nicolle PFS NEO1
# Figure 4b: LapNet PFS NEO1 (only SD/PR Pre samples)
# Figure 4c: LapNet OS NEO1 (only SD/PR Pre samples)

# Extended 2a: LapNet ESTIMATE
# Extended 2b: LapNet Post vs Pre, heatmap
# Extended 2c: LapNet Post vs Pre, boxplots of top DEGs
# Extended 2d: LapNet Post vs Pre (paired), GSEA dotplot + ridgeplot
# Extended 2e: Cassier Post vs Pre, core EMT
# Extended 2f: LapNet Post vs Pre, SD vs PR, core EMT

# Extended 4a: LapNet correlation NEO1:PFS
# Extended 4b: Nicolle correlation NEO1:PFS
# Extended 4c: Nicolle OS NEO1
# Extended 4d: LapNet PFS NEO1 (all samples Pre)
# Extended 4e: LapNet OS NEO1 (all samples Pre)


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
  library(hrbrthemes)
  library(ggpubr)
  library(msigdbr)
  library(patchwork)
  library(singscore)
  library(survival)
  library(survminer)
  library(tidyestimate)
})

source("R/helper.functions.R", echo=TRUE)

set.seed(1)

out.dir <- "results/"
data.dir <- "data/"


# Prepare data ----

## datasets ----

# Nicolle
load(file="data/LapNet/DGElist_GSE253260_LAPC_full.pdata.RData")
pdata_Nicolle <- y$samples
head(pdata_Nicolle)
rm(y)

# Linehan
load(file=paste0(data.dir, "DGElist_Linehan.Folfirinox_paired.RData"))
pdata_Linehan <- y$samples
pdata_Linehan$Pre.Post <- as.factor(pdata_Linehan$Pre.Post)
pdata_Linehan$Pre.Post <- relevel(pdata_Linehan$Pre.Post, ref = "Pre")
head(pdata_Linehan)
rm(y)

# LapNet
# using tumor and stromal samples for ESTIMATE
load(file=paste0(data.dir, "DGElist_all.RData"))
pdata_all <- y$samples
counts_all <- cpm(y, log=F, normalized.lib.sizes = T)
rownames(counts_all) <- y$genes$Symbol
rm(y)

# for all other plots, using only tumors
load(file=paste0(data.dir, "DGElist_norm.RData"))
pdata <- y$samples
#load(file=paste0(data.dir, "DGElist_tumor.all_full.metadata.RData"))
logcpm <- cpm(y, prior.count=2, log=TRUE, normalized.lib.sizes = T)
rownames(logcpm) <- y$genes$Symbol
top <- read.csv(file = paste0(out.dir, "DEGs_Post.vs.Pre.csv"))


## colors ----

condition.colors <- as.vector(alphabet())[c(25, 20)]
names(condition.colors) <- levels(pdata$timepoint)

response.colors <- as.vector(pals::alphabet())[c(1,2,7)]
pdata$Response <- as.factor(pdata$Response)
names(response.colors) <- levels(as.factor(pdata$Response))


## pdata subsets ----

pdata.pre <- pdata[pdata$timepoint=="Pre", ]
pdata.pre2 <- pdata.pre[pdata.pre$Response != "PD", ] # remove PD sample
paired_samples <- c("01-007","02-001","02-003","02-014")
pdata.paired <- pdata[pdata$ID_Lapnet %in% paired_samples, ]

pdata2 <- pdata[pdata$Response != "PD", ]
pdata2$Response <- factor(pdata2$Response)
pdata2$Response <- relevel(pdata2$Response, ref = "SD")


## NEO1 groups ----

# using median as cutoff
summary(pdata.pre$NEO1) # median = 6.245
pdata.pre$NEO1.group <- pdata.pre$NEO1 > quantile(pdata.pre$NEO1, 0.5)
pdata.pre$NEO1.group <- gsub("TRUE", "High", pdata.pre$NEO1.group)
pdata.pre$NEO1.group <- gsub("FALSE", "Low", pdata.pre$NEO1.group)
pdata.pre$NEO1.group <- as.factor(pdata.pre$NEO1.group)
pdata.pre$NEO1.group <- relevel(pdata.pre$NEO1.group, ref = "Low")
table(pdata.pre$NEO1.group)

pdata.pre2$NEO1.group <- pdata.pre2$NEO1 > quantile(pdata.pre2$NEO1, 0.5)
pdata.pre2$NEO1.group <- gsub("TRUE", "High", pdata.pre2$NEO1.group)
pdata.pre2$NEO1.group <- gsub("FALSE", "Low", pdata.pre2$NEO1.group)
pdata.pre2$NEO1.group <- as.factor(pdata.pre2$NEO1.group)
pdata.pre2$NEO1.group <- relevel(pdata.pre2$NEO1.group, ref = "Low")
table(pdata.pre2$NEO1.group)


# Figure 3 ----

## 3b ----

# select the maximum PValue corresponding to FDR < 0.05
pCutoff <- max(top$PValue[top$FDR < 0.05 & !is.na(top$FDR)])

EnhancedVolcano(top,
                lab = top$Symbol,
                x = 'logFC',
                y = 'PValue',
                pCutoff = pCutoff,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 3.0,
                max.overlaps = 20,
                legendPosition = "top",
                drawConnectors = T,
                boxedLabels = T,
                title = "Volcano plot",
                subtitle = bquote(italic("Post vs Pre"))
)

ggsave(paste0(out.dir, "Volcano_DEGs_Post.vs.Pre.png"), width = 8, height = 12)
ggsave(paste0(out.dir, "Volcano_DEGs_Post.vs.Pre.pdf"), width = 8, height = 12)


## 3c ----

em2 <- readRDS(paste0(out.dir, "GSEA_Hallmarks_Post.vs.Pre.rds"))

p1 <- dotplot(em2, x = "p.adjust", orderBy = "pvalue", showCategory = 10)
p2 <- ridgeplot(em2, orderBy = "pvalue", showCategory = 10) +
  theme(axis.text.y = element_blank()) +
  xlab("enrichment")

combined_plot <- (p1 + p2) + 
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = "GSEA for Hallmark Pathways",
    subtitle = "Post vs Pre",
    caption = "Data Source: [LapNET only Tumor tissues]"
  ) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 14))
combined_plot

ggsave(paste0(out.dir, "GSEA_Hallmarks_Post.vs.Pre.png"), combined_plot, width = 10, height = 8)
ggsave(paste0(out.dir, "GSEA_Hallmarks_Post.vs.Pre.pdf"), combined_plot, width = 10, height = 8)



## 3d ----

ggplot(data = pdata, aes(x = timepoint, y = core.EMT_score)) +
  geom_violin(aes(fill = timepoint), width = 1, alpha = 0.5, trim = F, position = position_dodge(width = 0.9)) +
  geom_boxplot(aes(fill = timepoint), width = 0.2, outlier.size = 0.1, outlier.colour = "black", alpha = 0.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = condition.colors) +
  stat_compare_means(method = "wilcox", label = "p.format", vjust = -4, 
                     comparisons = list(c("Pre","Post"))) +
  theme_minimal() +
  labs(y="Signature Score", x = "", title="Core EMT Signature") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))

ggsave(paste0(out.dir, "Core.EMT_Pre.vs.Post.png"), width = 6, height = 5)
ggsave(paste0(out.dir, "Core.EMT_Pre.vs.Post.pdf"), width = 6, height = 5)


## 3e ----

ggplot(data = pdata_Linehan, aes(x = Pre.Post, y = core.EMT_score)) +
  geom_violin(aes(fill = Pre.Post), width = 1, alpha = 0.5, trim = F, position = position_dodge(width = 0.9)) +
  geom_boxplot(aes(fill = Pre.Post), width = 0.2, outlier.size = 0.1, outlier.colour = "black", alpha = 0.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = condition.colors) +
  stat_compare_means(method = "wilcox", label = "p.format", vjust = -4, 
                     comparisons = list(c("Pre","Post"))) +
  theme_minimal() +
  labs(y="Signature Score", x = "", title="Core EMT Signature (Linehan-seq)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))

ggsave(paste0(out.dir, "Linehan_Core.EMT_Pre.vs.Post.png"), width = 6, height = 5)
ggsave(paste0(out.dir, "Linehan_Core.EMT_Pre.vs.Post.pdf"), width = 6, height = 5)


## 3f ----

em2 <- readRDS(paste0(out.dir, "GSEA_Hallmarks_Folfi.NP137.vs.Folfi_paired.interaction.rds"))

p1 <- dotplot(em2, x = "p.adjust", orderBy = "pvalue", showCategory = 10)
p2 <- ridgeplot(em2, orderBy = "pvalue", showCategory = 10) +
  theme(axis.text.y = element_blank()) +
  xlab("enrichment")

combined_plot <- (p1 + p2) + 
  plot_layout(guides = "collect") + # This collects all legends into one
  plot_annotation(
    title = "GSEA for Hallmark Pathways",
    subtitle = "FOLFIRINOX_NP137.vs.FOLFIRINOX (paired))",
    caption = "Data Source: [Integrated Lapnet-seq/Linehan-seq]"
  ) & 
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 14))
combined_plot

ggsave(paste0(out.dir, "GSEA_Hallmarks_Folfi.NP137.vs.Folfi_paired.interaction.png"), combined_plot, width = 10, height = 8)
ggsave(paste0(out.dir, "GSEA_Hallmarks_Folfi.NP137.vs.Folfi_paired.interaction.pdf"), combined_plot, width = 10, height = 8)


## 3g ----

ggplot(data = pdata.pre, aes(x = PFS_time, y = singscore_EMT)) +
  geom_point(color = "darkblue", alpha = 0.6) +
  geom_smooth(method = "lm", color = "red", fill = "pink", alpha = 0.3) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 5) +
  theme_minimal() +
  labs(x="PFS (days)", y = paste0("Signature Score"),
       title=paste0("full EMT Score vs PFS")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))

ggsave(paste0(out.dir, "full.EMT.vs.PFS.png"), width = 6, height = 5)
ggsave(paste0(out.dir, "full.EMT.vs.PFS.pdf"), width = 6, height = 5)


# Figure 4 ----

## 4a ----

# median PFS for NEO1 High is 194 days and for NEO1 Low is 259 days
quantile(pdata_Nicolle$exp_NEO1, 0.5) # 4.32116 
score.group <- pdata_Nicolle$exp_NEO1 > quantile(pdata_Nicolle$exp_NEO1, 0.5)
table(score.group)
fit.pfs <- survfit(Surv(PFS_time, PFS_event) ~ score.group, data = pdata_Nicolle)
p <- ggsurvplot(fit.pfs, pval = T, risk.table = T,
           palette = c("darkorange","darkorchid"),
           fontsize = 4,  font.main = c(16, "bold", "black"),
           surv.median.line = "hv",
           break.time.by = 200,
           pval.size = 4, pval.method = T,
           legend = "top", 
           legend.title = "expression", 
           legend.labs = c("Low", "High"),
           title = paste0("PFS by NEO1", " expression"),
           subtitle = "Dataset GSE253260"
)
ggsave(paste0(out.dir, "PFS_NEO1.group_Nicolle.pdf"), p$plot, width = 7, height = 6)



## 4b ----
# no PD sample

fit.pfs <- survfit(Surv(PFS_time, PFS_status) ~ NEO1.group, data = pdata.pre2)
#png(paste0(out.dir, "PFS_NEO1.group_no.PD.png"), width = 600, height = 600)
p <- ggsurvplot(fit.pfs, pval = T, risk.table = T,
                palette = c("cornflowerblue","brown1"),
                fontsize = 4,  font.main = c(16, "bold", "black"),
                surv.median.line = "hv",
                break.time.by = 100,
                pval.size = 5, pval.method = T,
                legend = "top", 
                legend.title = "score", 
                legend.labs = c("Low", "High"),
                title = paste0("PFS by NEO1 Group"),
                subtitle = "NEO1 expression group in tumor samples Pre-treatment"
)
ggsave(paste0(out.dir, "PFS_NEO1.group_no.PD.pdf"), p$plot, width = 7, height = 6)

# median survival for NEO1 High is 476 days and for NEO1 Low is 311 days


## 4c ----

fit.os <- survfit(Surv(OS_time, OS_status) ~ NEO1.group, data = pdata.pre2)
#png(paste0(out.dir, "OS_NEO1.group_no.PD.png"), width = 600, height = 600)
p <- ggsurvplot(fit.os, pval = T, risk.table = T,
                palette = c("cornflowerblue","brown1"),
                fontsize = 4,  font.main = c(16, "bold", "black"),
                surv.median.line = "hv",
                break.time.by = 100,
                pval.size = 5, pval.method = T,
                legend = "top", 
                legend.title = "score", 
                legend.labs = c("Low", "High"),
                title = paste0("OS by NEO1 Group"),
                subtitle = "NEO1 expression group in tumor samples Pre-treatment"
)
ggsave(paste0(out.dir, "OS_NEO1.group_no.PD.pdf"), p$plot, width = 7, height = 6)

# median survival for NEO1 High is NA days and for NEO1 Low is 501 days



# Extended 2 ----

## E2a ----

# using tumor and non-tumor samples
filtered <- filter_common_genes(counts_all, 
                                id = "hgnc_symbol", 
                                tidy = FALSE, 
                                tell_missing = TRUE, 
                                find_alias = TRUE)
scored <- estimate_score(filtered,
                         is_affymetrix = TRUE)
head(scored)
plot_purity(scored, is_affymetrix = TRUE)

ggsave(paste0(out.dir, "LapNet_ESTIMATE_purity_all.samples.png"), width=10, height=10)
ggsave(paste0(out.dir, "LapNet_ESTIMATE_purity_all.samples.pdf"), width=10, height=10)

pdata.tumor <- pdata_all[pdata_all$Structures=="Tumor", ]
score.tumor <- scored[scored$sample %in% rownames(pdata.tumor), ]
all(rownames(pdata.tumor)==score.tumor$sample) # TRUE
score.tumor$sample <- pdata.tumor$ID_Lapnet
options(ggrepel.max.overlaps = Inf)
plot_purity(score.tumor, is_affymetrix = TRUE)

ggsave(paste0(out.dir, "LapNet_ESTIMATE_purity_tumor.samples.png"), width=8, height=8)
ggsave(paste0(out.dir, "LapNet_ESTIMATE_purity_tumor.samples.pdf"), width=8, height=8)


## E2b ----

top50 <- top[top$FDR < 0.01 & top$logFC>abs(1), ]
Sel <- logcpm[rownames(logcpm)%in%top50$Symbol, ]

#png(filename = paste0(out.dir, "heatmap_Post.vs.Pre_FDR.001_FC.2.png"), width = 600, height = 900)
pdf(file = paste0(out.dir, "heatmap_Post.vs.Pre_FDR.001_FC.2.pdf"), width = 5, height = 8)
aheatmap(Sel, scale = "row", distfun = "euclidean",
         annCol = list(timepoint = pdata$timepoint), 
         annColors = list(condition.colors),
         annLegend = T,
         main="Post vs Pre\n[FDR<0.01, FC>|2|)")
dev.off()



## E2c ----

# select the top 12 genes with lowest and highest fold change
sel.genes <- c(rev(na.omit(top$Symbol[top$logFC < 0])[1:6]),
               na.omit(top$Symbol[top$logFC > 0])[1:6])
all(rownames(pdata)==colnames(logcpm))
plot_data <- cbind(pdata, t(logcpm[sel.genes, ]))
plot_data_facet <- plot_data %>%
  select(timepoint, ID_Lapnet, all_of(sel.genes)) %>%
  pivot_longer(cols = all_of(sel.genes),
               names_to = "gene",
               values_to = "expression") %>%
  mutate(gene = factor(gene, levels = sel.genes)) # Convert 'gene' to ordered factor

boxplot_facet <- ggplot(plot_data_facet, aes(x = timepoint, y = expression, fill = timepoint)) +
  geom_boxplot() +
  scale_fill_manual(values = condition.colors) +
  labs(y = "expression [logcpm]", fill = "Condition") +
  theme_bw() +
  xlab("") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 14), # Increase the size here (e.g., to 14)
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ gene, ncol = 4)
print(boxplot_facet)

ggsave(paste0(out.dir, "boxplot_facet_DEGs_Post.vs.Pre.png"), plot = boxplot_facet, width = 7, height = 7)
ggsave(paste0(out.dir, "boxplot_facet_DEGs_Post.vs.Pre.pdf"), plot = boxplot_facet, width = 7, height = 7)


## E2d ----

em2 <- readRDS(paste0(out.dir, "GSEA_Hallmarks_Post.vs.Pre_paired.rds"))

p1 <- dotplot(em2, x = "p.adjust", orderBy = "pvalue", showCategory = 10)
p2 <- ridgeplot(em2, orderBy = "pvalue", showCategory = 10) +
  theme(axis.text.y = element_blank()) +
  xlab("enrichment")

combined_plot <- (p1 + p2) + 
  plot_layout(guides = "collect") + # This collects all legends into one
  plot_annotation(
    title = "GSEA for Hallmark Pathways",
    subtitle = "Post vs Pre (paired)",
    caption = "Data Source: [LapNET only Tumor tissues]"
  ) & # The '&' operator applies the theme/setting to all subplots
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 14))
combined_plot

ggsave(paste0(out.dir, "GSEA_Hallmarks_Folfi.NP137.vs.Folfi_paired.interaction.png"), combined_plot, width = 10, height = 8)


## E2e ----

library(Seurat)

# data from Cassier et al.
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225691
load(paste0(data.dir, "Cassier.tumor.RData"))

colors <- as.vector(alphabet(n=24))
cluster.colors = colors[1:12]
condition.colors.sc = c("brown","orange")
DimPlot(tumor)
head(tumor)
table(tumor$ID)
tumor$ID <- as.factor(tumor$ID)

gsea.results <- read.csv(file = paste0(out.dir, "/GSEA_Hallmarks_Post.vs.Pre_v2.csv"), row.names = 1)
head(gsea.results)
core.genes <- gsea.results["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "core_enrichment"]
core.genes <- unlist(strsplit(core.genes, "/"))
sort(core.genes) # 45 genes

intersect(core.genes, rownames(tumor))
setdiff(core.genes, rownames(tumor)) # MATN3

tumor <- AddModuleScore(tumor, 
                        name = "Core.EMT_Score",
                        list(core.genes))

p <- VlnPlot(tumor, features = "Core.EMT_Score1", 
             group.by = "ID", cols = condition.colors.sc) +
  labs(title = "Core EMT Score",
       subtitle = "Cassier et al. (tumor cells)") + xlab("") +
  stat_compare_means(comparisons = list(c("C1D1","C3D1")), label = "p.format", label.y.npc = "top") +
  theme(axis.text.x = element_text(angle = 0, size = 16, hjust=0.5)) +
  ylim(-0.3, 0.7)

ggsave(paste0(out.dir, "Cassier_Core.EMT_Score.png"), p, width = 10, height = 5)
ggsave(paste0(out.dir, "Cassier_Core.EMT_Score.pdf"), p, width = 10, height = 5)



## E2f ----

ggplot(data = pdata2, aes(x = Response, y = core.EMT_score)) +
  geom_violin(aes(fill = timepoint), width = 1, alpha = 0.5, trim = F, 
              scale = "width", position = position_dodge(width = 1)) +
  geom_boxplot(aes(fill = timepoint), width = 0.2, outlier.size = 0.1, outlier.colour = "black", 
               alpha = 0.2, position = position_dodge(width = 1)) +
  scale_fill_manual(values = condition.colors) +
  stat_compare_means(method = "wilcox", label = "p.format", vjust = -4, aes(group = timepoint)) +
  theme_minimal() +
  labs(y="Signature Score", x = "", title="Core EMT Signature") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))

ggsave(paste0(out.dir, "Core_EMT_Pre.vs.Post_Response.png"), width = 8, height = 5)
ggsave(paste0(out.dir, "Core_EMT_Pre.vs.Post_Response.pdf"), width = 8, height = 5)


# Extended 4 ----

## E4a ----

ggplot(data = pdata.pre2, aes(x = PFS_time, y = NEO1)) +
  geom_point(color = "darkblue", alpha = 0.6) +
  geom_smooth(method = "lm", color = "red", fill = "pink", alpha = 0.3) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 5) +
  theme_minimal() +
  labs(x="PFS (days)", y = paste0("NEO1 expression (logCPM)"),
       title=paste0("NEO1 expression vs PFS")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))

ggsave(paste0(out.dir, "NEO1.vs.PFS_no.PD.png"), width = 6, height = 5)
ggsave(paste0(out.dir, "NEO1.vs.PFS_no.PD.pdf"), width = 6, height = 5)


## E4b ----

p <- ggplot(data = pdata_Nicolle, aes(x = PFS_time, y = exp_NEO1)) +
  geom_point(color = "steelblue", alpha = 0.6) +
  geom_smooth(method = "lm", color = "mediumseagreen", fill = "lightgreen", alpha = 0.3) +
  stat_cor(method = "pearson", label.x.npc = "center", label.y.npc = "bottom", size = 5) +
  theme_minimal() +
  labs(x="Progression Free Survival (days)", y = paste0("NEO1 expression (logCPM)"),
       title=paste0("mFOLFIRINOX"),
       subtitle = "") + # "Dataset GSE253260") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0, size = 16, face = "bold"))

ggsave(paste0(out.dir, "Nicolle_PFS_NEO1_expression.png"), p, width=6, height=5)
ggsave(paste0(out.dir, "Nicolle_PFS_NEO1_expression.pdf"), p, width=6, height=5)


## E4c ----

# median OS for NEO1 High is 443 days and for NEO1 Low is 438 days
quantile(pdata_Nicolle$exp_NEO1, 0.5) # 4.32116 
score.group <- pdata_Nicolle$exp_NEO1 > quantile(pdata_Nicolle$exp_NEO1, 0.5)
table(score.group)

fit.os <- survfit(Surv(OS_time, OS_event) ~ score.group, data = pdata_Nicolle)
ggsurvplot(fit.os, pval = T, risk.table = T,
           palette = c("darkorange","darkorchid"),
           fontsize = 4,  font.main = c(16, "bold", "black"),
           surv.median.line = "hv",
           break.time.by = 200,
           pval.size = 4, pval.method = T,
           legend = "top", 
           legend.title = "expression", 
           legend.labs = c("Low", "High"),
           title = paste0("OS by NEO1", " expression"),
           subtitle = "Dataset GSE253260"
)
ggsave(paste0(out.dir, "OS_NEO1.group_Nicolle.pdf"), p$plot, width = 7, height = 6)


## E4d ----

fit.pfs <- survfit(Surv(PFS_time, PFS_status) ~ NEO1.group, data = pdata.pre)
#png(paste0(out.dir, "PFS_NEO1.group.png"), width = 600, height = 600)
p <- ggsurvplot(fit.pfs, pval = T, risk.table = T,
           palette = c("cornflowerblue","brown1"),
           fontsize = 4,  font.main = c(16, "bold", "black"),
           surv.median.line = "hv",
           break.time.by = 100,
           pval.size = 5, pval.method = T,
           legend = "top", 
           legend.title = "score", 
           legend.labs = c("Low", "High"),
           title = paste0("PFS by NEO1 Group"),
           subtitle = "NEO1 expression group in tumor samples Pre-treatment"
)

ggsave(paste0(out.dir, "PFS_NEO1.group.pdf"), p$plot, width = 7, height = 6)

# median survival for NEO1 High is 476 days and for NEO1 Low is 311 days


## E4e ----

fit.os <- survfit(Surv(OS_time, OS_status) ~ NEO1.group, data = pdata.pre)
#png(paste0(out.dir, "OS_NEO1.group.png"), width = 600, height = 600)
p <- ggsurvplot(fit.os, pval = T, risk.table = T,
           palette = c("cornflowerblue","brown1"),
           fontsize = 4,  font.main = c(16, "bold", "black"),
           surv.median.line = "hv",
           break.time.by = 100,
           pval.size = 5, pval.method = T,
           legend = "top", 
           legend.title = "score", 
           legend.labs = c("Low", "High"),
           title = paste0("OS by NEO1 Group"),
           subtitle = "NEO1 expression group in tumor samples Pre-treatment"
)

ggsave(paste0(out.dir, "OS_NEO1.group.pdf"), p$plot, width = 7, height = 6)

# median survival for NEO1 High is NA days and for NEO1 Low is 501 days




# end ----
sessionInfo()

