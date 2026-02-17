
# main --------------------------------------------------------------------

# microBulk RNAseq: human FFPE + LCM, PDAC tissues
# paired fastq files
# checked with fastQC + multiqc
# trimmed with Trim-Galore!

# Initial code based on:
# https://bioconductor.org/packages/release/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html#read-alignment-and-quantification

# considerations for microbulk + FFPE:
# expect high level of duplicates and lower quality in fastQC
# stringent trimming may be necessary
# more strict feature counting
# expect low mapping and assignment rates
# loss of reads may be compensated by high sequencing depth
# favor pathway analyses over single-gene comparisons


# Settings ---------------------------------------------------------------

rm(list=ls())
gc()

suppressPackageStartupMessages({
  library(Rsubread)
  library(Rsamtools)
  library(stringr)
  library(edgeR)
  library(org.Hs.eg.db)
  library(ggfortify)
  library(eulerr)
  library(pals)
  library(gridExtra)
  library(EnhancedVolcano)
  library(pals)
  library(scales)
  library(dplyr)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(hrbrthemes)
  library(ggpubr)
  library(patchwork)    
})

source("R/helper.functions.R", echo=TRUE)

set.seed(1)

out.dir <- "results/"
data.dir <- "data/"


# loading + inspection + alignment ----------------------------------------

fastq.files <- list.files(path = "trimmed_fastq/", pattern = ".gz", recursive = F, full.names = F)
full.path <- list.files(path = "trimmed_fastq/", pattern = ".gz", recursive = F, full.names = T)
samples <- str_sub(basename(fastq.files), 1, -20)

# Build a genome index
# buildindex(basename = "hg38", reference = "hg38.fa")

# Align reads
forward <- full.path[c(TRUE, FALSE)]
reverse <- full.path[c(FALSE, TRUE)]
all.bam <- make.names(paste0(unique(samples), ".bam"))
align(index="hg38", readfile1= forward, readfile2= reverse,
      type = "rna", nthreads = 12, annot.inbuilt = "hg38",
      input_format="gzFASTQ", output_file= all.bam)

# moved bams to separate directory
all.bam <- list.files(path = "bam/", pattern = "\\.bam$", recursive = T, full.names = T)

# summarize the files created during aligment:
file_list <- list.files(path = "bam/", pattern = "\\.summary$", full.names = TRUE)
# Read the first file to initialize the data frame
prop.mapped.table <- read.table(file_list[1], sep = "\t", header = FALSE, stringsAsFactors = FALSE)
prop.mapped.table <- as.data.frame(t(prop.mapped.table))
colnames(prop.mapped.table) <- prop.mapped.table[1, ]
prop.mapped.table <- prop.mapped.table[-1, ]
# Add a column for the file name
prop.mapped.table$File <- basename(file_list[1])
# Loop through the rest of the files and append them to the data frame
for (i in 2:length(file_list)) {
  temp_data <- read.table(file_list[i], sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  temp_data <- as.data.frame(t(temp_data))
  colnames(temp_data) <- temp_data[1, ]
  temp_data <- temp_data[-1, ]
  temp_data$File <- basename(file_list[i])
  prop.mapped.table <- rbind(prop.mapped.table, temp_data)
}
head(prop.mapped.table)

write.csv(prop.mapped.table, file = "prop.mapped.csv")

# Quantify read counts
fc <- featureCounts(all.bam, annot.inbuilt="hg38", 
                    nthreads= 12, 
                    isPairedEnd=T, 
                    useMetaFeatures = TRUE,
                    countMultiMappingReads = F,
                    strandSpecific = 0, # checked options 1 and 2 with even less alignments
                    minMQS = 20, # default parameter is zero
                    requireBothEndsMapped = T # increase stringency for FFPE samples
)
head(fc$counts)

save(fc, file= paste0(data.dir, "counts.RData"))
write.csv(prop.mapped.table, file = paste0(data.dir, "prop.mapped.table.csv"))
# save raw counts for GEO
write.csv(fc$counts, file = paste0(data.dir, "raw_counts.csv"))


## QC -----------------------------------------------------

load(paste0(data.dir, "counts.RData"))
head(fc$counts)

# check stats

names(fc)
head(fc$stat)
stat <- fc$stat
rownames(stat) <- stat$Status
stat <- stat[,-1]
head(stat)

prop.mapped <- read.csv("data/LapNet/prop.mapped.csv", row.names=NULL)
prop.mapped$NumMapped <- prop.mapped$Mapped_fragments
prop.mapped$NumTotal <- prop.mapped$Total_no.trimming
head(prop.mapped)

all(prop.mapped$bam.file==colnames(stat))
prop.mapped <- cbind(prop.mapped, t(stat))
prop.mapped$prop <- prop.mapped$NumMapped/prop.mapped$NumTotal
prop.mapped$prop.assigned <- prop.mapped$Assigned/prop.mapped$NumTotal
prop.mapped$prop.untrimmed <- prop.mapped$Total_fragments/prop.mapped$NumTotal
rownames(prop.mapped) <- str_sub(rownames(prop.mapped), 1, -5)
head(prop.mapped)

p1 <- ggplot(prop.mapped) + 
  geom_bar(aes(fill="Total reads", y=NumTotal, x=sample), position="stack", stat="identity") +
  geom_bar(aes(fill="Reads untrimmed", y=Total_fragments, x=sample), position="stack", stat="identity") +
  geom_bar(aes(fill="Reads mapped", y=NumMapped, x=sample), position="stack", stat="identity") +
  geom_bar(aes(fill="Reads assigned", y=Assigned, x=sample), position="stack", stat="identity") +
  #  geom_point(aes(fill="Duplication", y=av_dups, x=rownames(prop.mapped))) + 
  scale_fill_brewer(palette = "Blues", direction = -1) +
  xlab("") + ylab("Number of reads") +
  ggtitle("Mapping QC - Number") +
  theme_minimal() + labs(fill = "") +
  theme(axis.text.x = element_text(size = 12,angle = 90, vjust = 0.5, hjust=1),
        title = element_text(size=14, face="bold"),
        legend.text = element_text("", size=12))

p2 <- ggplot(prop.mapped) + 
  geom_bar(aes(fill="Total reads", y=1, x=sample), position="stack", stat="identity") +
  geom_bar(aes(fill="Reads untrimmed", y=prop.untrimmed, x=sample), position="stack", stat="identity") +
  geom_bar(aes(fill="Reads mapped", y=prop, x=sample), position="stack", stat="identity") +
  geom_bar(aes(fill="Reads assigned", y=prop.assigned, x=sample), position="stack", stat="identity") +
  #  geom_point(aes(fill="Duplication", y=av_dups, x=rownames(prop.mapped))) + 
  scale_fill_brewer(palette = "Blues", direction = -1) +
  xlab("") + ylab("Proportion of reads") +
  ggtitle("Mapping QC - Proportion") +
  theme_minimal() + labs(fill = "") +
  theme(axis.text.x = element_text(size = 12,angle = 90, vjust = 0.5, hjust=1),
        title = element_text(size=14, face="bold"),
        legend.text = element_text("", size=12))

p1/p2
ggsave(paste0(out.dir, "mapping_QC_barplots.png"),
       width = 15, height = 15)



# edgeR preprocessing -----------------------------------------------------

load(paste0(data.dir, "counts.RData"))
head(fc$counts)

pdata <- read.csv(paste0(data.dir, "pdata_LapNet.csv"))
head(pdata)
table(pdata$timepoint)
# Post  Pre 
#  9   27 
table(pdata$Structures)
# Stroma  Tumor 
#   6     30 

rownames(pdata) <- pdata$bam.file
pdata <- pdata[colnames(fc$counts), ]

all(colnames(fc$counts) == pdata$bam.file)

# change names
rownames(pdata) <- pdata$R_name
colnames(fc$counts) <- rownames(pdata)

y <- DGEList(counts= fc$counts, samples = pdata)
dim(y) # 28395    36

# aggregate counts of technical replicate
y <- sumTechReps(y, ID = y$samples$ID)
dim(y) # 28395    35

# Add gene annotation
y$genes <- fc$annotation[, c("Length","GeneID"), drop=FALSE]
y$genes$Symbol <- mapIds(org.Hs.eg.db, rownames(y), keytype="ENTREZID", column="SYMBOL")
y$genes$EntrezID <- rownames(y)
head(y$genes)
y <- y[!is.na(y$genes$Symbol), ]
dim(y) # 27773    35

save(y, file= paste0(data.dir, "DGElist_all.RData"))


# inspection -----------------------------------------------------

load(file= paste0(data.dir, "DGElist_all.RData"))
pdata <- y$samples
head(pdata)

# extract counts, cpm and logcpm for different metrics
counts <- y$counts
rownames(counts) <- y$genes$Symbol
counts[1:5,1:5]
dim(counts)
logcpm <- cpm(y, prior.count=2, log=T, normalized.lib.sizes = T)

par(mar=c(10,5,5,5), mfrow=c(3,1))

boxplot(logcpm, las = 2, ylab = "logcpm") # 
barplot(colSums(counts), las = 2, main = "Read counts per sample")
barplot(rowSums(counts), las = 2, main = "Read counts per feature")
# bad QC sample X24N0_2454.06.006 has 1072 reads
# others range between 132548 and 12813357 reads 
# about half of them have at least 1M reads

sort(rowSums(counts), decreasing = T)[1:100]
sort(colSums(counts), decreasing = T)

colnames(counts)
sort(counts[, 1], decreasing = T)[1:20]
sort(counts[, 2], decreasing = T)[1:20]
sort(counts[, 3], decreasing = T)[1:20]
sort(counts[, 4], decreasing = T)[1:20]
sort(counts[, 5], decreasing = T)[1:20]
sort(counts[, 6], decreasing = T)[1:20]
sort(counts[, 7], decreasing = T)[1:20]

# usual "flagged" suspects, including MALAT1 and MT- genes
# not filtering as this seems similar across samples

show_col(alphabet())
condition.colors <- as.vector(pals::alphabet())[c(25,19)]
names(condition.colors) <- levels(pdata$timepoint)

logcpm.pre <- cpm(y, prior.count=2, log=TRUE, normalized.lib.sizes = F)
logcpm.pre[1:5,1:5]

par(mar=c(15,5,5,5), mfrow=c(1,1))
col.group <- as.factor(pdata$timepoint)
boxplot(logcpm.pre, las = 2, cex.axis = 1,
        main= "Data distribution before normalization", 
        ylab = "logcpm", col = condition.colors[col.group])


# PCA pre-normalization

t.counts <- t(logcpm.pre)
t.counts[1:5,1:10]

pca_res <- prcomp(t.counts, scale. = F)

p1 <- autoplot(pca_res, 
               data = pdata, 
               colour = 'timepoint', 
               size = 4, frame = F,
               label=F, label.size = NA, frame.type = 'norm') + #, label=F, label.size = NA, frame = F, shape = T) +
  theme_minimal() +
  scale_color_manual(values = condition.colors) +
  ggtitle("PCA before normalization",
          subtitle = "by Condition")
p2 <- autoplot(pca_res, 
               data = pdata, 
               colour = 'Structures', 
#               fill = 'Structures', 
size = 4, frame = F,
label=F, label.size = NA, frame.type = 'norm') + #, label=F, label.size = NA, frame = F, shape = T) +
  theme_minimal() +
  scale_color_manual(values = c("orange","navyblue")) +
  ggtitle("PCA before normalization",
          subtitle = "by Tissue")
p3 <- autoplot(pca_res, 
               data = pdata, 
               colour = 'ID_Patho', 
               size = 4, frame = F,
               label=F, label.size = NA, frame.type = 'norm') + #, label=F, label.size = NA, frame = F, shape = T) +
  theme_minimal() +
  scale_color_manual(values = as.vector(scPalette(length(unique(pdata$ID_Patho))))) +
  guides(colour = "none") +
  ggtitle("PCA before normalization",
          subtitle = "by Patient ID")

grid.arrange(p1, p2, p3, ncol=3)



## normalization ----

# Model
design <- model.matrix(~1, data = pdata)
ncol(design)
head(design)

# Filtering
keep <- filterByExpr(y, design)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y) # 

# Normalization for composition bias
y <- calcNormFactors(y)

logcpm <- cpm(y, prior.count=2, log=TRUE, normalized.lib.sizes = T)
logcpm[1:5,1:5]

par(mar=c(15,5,5,5), mfrow=c(1,1))
col.group <- as.factor(pdata$timepoint)
boxplot(logcpm.pre, las = 2, cex.axis = 1,
        main= "Data distribution after normalization", 
        ylab = "logcpm", col = condition.colors[col.group])


# PCA post-normalization

t.counts <- t(logcpm)
t.counts[1:5,1:10]
pca_res2 <- prcomp(t.counts, scale. = F)

p1 <- autoplot(pca_res2, 
               data = pdata, 
               colour = 'timepoint', 
               size = 4, frame = F,
               label=F, label.size = NA, frame.type = 'norm') + #, label=F, label.size = NA, frame = F, shape = T) +
  theme_minimal() +
  scale_color_manual(values = condition.colors) +
  ggtitle("PCA after normalization",
          subtitle = "by Condition")
p2 <- autoplot(pca_res2, 
               data = pdata, 
               colour = 'Structures', 
               #               fill = 'Structures', 
               size = 4, frame = F,
               label=F, label.size = NA, frame.type = 'norm') + #, label=F, label.size = NA, frame = F, shape = T) +
  theme_minimal() +
  scale_color_manual(values = c("orange","navyblue")) +
  ggtitle("PCA after normalization",
          subtitle = "by Tissue")
p3 <- autoplot(pca_res2, 
               data = pdata, 
               colour = 'ID_Patho', 
               size = 4, frame = F,
               label=F, label.size = NA, frame.type = 'norm') + #, label=F, label.size = NA, frame = F, shape = T) +
  theme_minimal() +
  scale_color_manual(values = as.vector(scPalette(length(unique(pdata$ID_Patho))))) +
  guides(colour = "none") +
  ggtitle("PCA after normalization",
          subtitle = "by Patient ID")

grid.arrange(p1, p2, p3, ncol=3)



# end ---------------------------------------------------------------------
sessionInfo()





