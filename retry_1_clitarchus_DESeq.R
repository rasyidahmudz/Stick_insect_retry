### DEseq Script ###
# Libraries
library(Rsubread)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(Rgraphviz)
library(dplyr)
library(tidyr)

# IMPT FILES
setwd("~/Desktop/Stick_insect")
bam_dir <- "/Users/rm2344/Desktop/Stick_insect/clitarchus_bamfiles"
bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)
annotation <- "/Users/rm2344/Desktop/Stick_insect/cli_annevo_gemoma_combined.gtf"

### PATTERNS OF GENE EXPRESSION IN SAMPLES ###
## GENERATE GENE COUNTS FROM BAM FILES ##
fc <- featureCounts(
  files = bam_files,
  annot.ext = annotation,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id", # before this was GTF.attrType = "Parent" for gff3
  useMetaFeatures = TRUE,
  isPairedEnd = TRUE,
  nthreads = 10
)
# shortcut from saved
fc <- readRDS("~/Desktop/Stick_insect/rds/featureCounts_fc.rds")

counts <- fc$counts # extract count matrix
colnames(counts) <- gsub("_Aligned.sortedByCoord.out.bam", "", colnames(counts)) # clean column names
# Pad numbers with zeros E.g. CLH_1L --> CLH_001L
nums <- as.numeric(sub("CLH_(\\d+)L", "\\1", colnames(counts)))
colnames(counts) <- sprintf("CLH_%03dL", nums)
# Reorder columns alphabetically (which now matches numeric order)
counts <- counts[, order(colnames(counts))]

# Create sample colData
sampledata <- data.frame(
  acclimation.temp = factor(c(
    20,20,20,20,20,20,20,20,20,20,20,20,
    10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
    20,20,20,20,20,20,20,
    10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
    20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,
    10,10,10,
    20
  )),
  test.temp = factor(c(
    32,32,20,20,32,15,15,5,5,5,
    20,15,32,32,32,20,15,15,5,5,
    5,10,10,25,25,25,30,30,30,30,
    25,25,25,30,25,25,5,5,15,15,
    20,20,32,10,10,10,25,30,25,32,
    30,35,35,35,35,35,32,5,5,15,
    15,20,32,10,10,10,10,35,20,30,
    30,35,35,10,35,35,15,20,20,30
  )),
  mass = c(
    1.347, 1.195, 0.998, 1.144, 1.001, 1.134, 1.18, 1.155, 1.348, 1.075,
    0.657, 0.981, 0.816, 1.432, 1.162, 1.328, 1.175, 1.155, 1.359, 1.327,
    1.087, 1.334, 0.905, 0.79, 1.233, 1.768, 1.258, 1.666, 1.421, 1.261,
    1.208, 1.088, 0.999, 1.194, 0.882, 1.145, 0.765, 1.465, 1.505, 1.614,
    0.785, 0.86, 1.52, 1.028, 1.155, 1.131, 1.441, 1.337, 1.35, 1.231,
    1.419, 1.26, 1.37, 1.373, 1.297, 1.13, 0.634, 1.489, 1.025, 1.185,
    0.983, 0.985, 1.216, 1.141, 0.686, 1.184, 1.129, 1.167, 1.159, 0.98,
    1.072, 1.071, 1.364, 1.218, 1.237, 1.471, 0.906, 0.693, 1.078, 1.006
  ),
  leaf.number = factor(c(
    1,2,3,4,5,6,7,8,9,10,
    11,12,14,15,16,17,18,20,21,22,
    23,24,25,26,27,28,29,30,31,36,
    37,39,40,48,45,46,
    1,2,3,4,5,6,7,8,35,36,
    37,38,39,40,42,43,44,45,46,47,
    9,10,11,12,13,14,15,49,50,51,
    52,53,54,55,56,57,58,59,60,61,
    1,3,5,7
  )),
  leaf.area.diff = c(
    0,3.195,32.3,57.217,41.829,55.95,45.989,27.799,1.138,26.011,
    33.059,37.65,0,40.825,0,75.39,23.235,33.197,0,0,
    0,43.257,34.221,31.281,30.513,72.321,14.793,27.348,58.911,78.167,
    71.752,58.822,59.038,45.504,48.148,35.985,
    19.345,15.789,26.582,51.469,50.859,37.747,18.926,14.679,37.13,15.817,
    70.699,42.026,104.305,0,46.162,0,0,0,0,0,
    29.413,0,2.281,34.275,47.817,50.562,78.414,36.379,23.1,24.678,
    29.614,0,69.337,65.262,44.088,0,0,27.146,0,0,
    38.979,44.857,20.068,44.746
  )
)
# Assign rownames AFTER creating metadata
rownames(sampledata) <- colnames(counts)

# Set as factors
sampledata$acclimation.temp <- factor(sampledata$acclimation.temp)
sampledata$test.temp <- factor(sampledata$test.temp)
sampledata$mass <- as.numeric(sampledata$mass)

# Sanity checks!
all(colnames(counts) == rownames(sampledata)) # check alignment

## CREATE DESeq2 OBJECT ##
# Interaction design
library(DESeq2)
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sampledata,
  design = ~ acclimation.temp * test.temp + mass
)
# Take a look:
dds

## NORMALIZATION ##

# Filtering Step #
# Determine whether genes meet our threshold
keep <- rowSums(counts(dds) >= 5 ) >= 3 # Check if 5 reads are too low for high depth reads

# Filter the dds object to keep genes that meet the threshold
dds_filt <- dds[keep,]

# Run DEseq on filtered dataset
library(DESeq2)
dds_filt <- DESeq(dds_filt)

# Check what contrasts are even possible
resultsNames(dds_filt)

# Normalize counts within filtered DDS object #
# Estimate the size factors
dds_filt <- estimateSizeFactors(dds_filt)
sizeFactors(dds_filt) # Check the size factors 

# Subset by ACCLIMATION TEMP
dds_10 <- dds_filt[, dds_filt$acclimation.temp == "10"]
dds_20 <- dds_filt[, dds_filt$acclimation.temp == "20"]

## UNSUPERVISED CLUSTERING ##
# Transform the data using the variance-stabilizing transformation
vst <- vst(dds_filt, blind=TRUE)

### PCA ###
# Extract PCA data
pcaData <- plotPCA(vst, intgroup = c("acclimation.temp", "test.temp", "mass"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# [ALL] plotPCA with ggplot2 by ACCLIMATION TEMP
library(ggplot2)
ggplot(pcaData, aes(PC1, PC2, color = acclimation.temp)) +
  geom_point(size = 3, alpha = 0.9) +
  
  scale_color_brewer(
    palette = "Dark2",
    name = "Acclimation temperature",
    labels = c(
      "10" = "10 degrees",
      "20" = "20 degrees"
    )
  ) +
  
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of samples by Acclimation temperature") +
  
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# [ALL] plotPCA with ggplot2 by TEST TEMP
ggplot(pcaData, aes(PC1, PC2, color = test.temp)) +
  geom_point(size = 3, alpha = 0.9) +
  
  scale_color_brewer(
    palette = "Dark2",
    name = "Test temperature",
    labels = c(
      "5" = "5 deg",
      "10" = "10 deg",
      "15" = "15 deg",
      "20" = "20 deg",
      "25" = "25 deg",
      "30" = "30 deg",
      "32" = "32 deg",
      "35" = "35 deg"
    )
  ) +
  
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of samples by Test temperature") +
  
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# [ALL] plotPCA with ggplot2 by MASS
ggplot(pcaData, aes(PC1, PC2, color = mass)) +
  geom_point(size = 3) +
  scale_color_viridis_c() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of samples by Mass") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# [ALL] plotPCA with ggplot2 by BATCH
ggplot(pcaData, aes(PC1, PC2, color = batch)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_brewer(
    palette = "Dark2",
    name = "Batch",
    labels = c(
      "1" = "1",
      "2" = "2",
      "3" = "3",
      "4" = "4",
      "0" = "NA"
    )
  ) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of samples by Batch") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# Multi-factor interpretation
pcaData$group <- interaction(
  pcaData$acclimation.temp,
  pcaData$test.temp
)
ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_point(size = 3, alpha = 0.9) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of samples by acclimation.temp_test.temp") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# Combine indiv plots into one figure
library(patchwork)
(plot_S1 | plot_S2) /
  (plot_S3 | plot_A)

### DIFFERENTIAL EXPRESSION ANALYSIS ###
# Note: Ref is acclimation.temp=10 and test.temp=5

# Interaction heatmap -plasticity landscape
int_terms <- c(
  "acclimation.temp20.test.temp10",
  "acclimation.temp20.test.temp15",
  "acclimation.temp20.test.temp20",
  "acclimation.temp20.test.temp25",
  "acclimation.temp20.test.temp30",
  "acclimation.temp20.test.temp32",
  "acclimation.temp20.test.temp35"
)
# Extract all
res_list <- lapply(int_terms, function(nm){
  results(dds_filt, name = nm)
})
# Build log2FC matrix
lfc_mat <- sapply(res_list, function(res) res$log2FoldChange)
colnames(lfc_mat) <- c("10", "15", "20", "25", "30", "32", "35")
# Select most informative genes
top_genes <- order(rowMeans(abs(lfc_mat), na.rm=TRUE), decreasing=TRUE)[1:200]
mat <- lfc_mat[top_genes, ]
# z-score per gene
mat_z <- t(scale(t(mat)))
# Heatmap
library(pheatmap)
pheatmap(
  mat_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  main = "Interaction effects: acclimation temp × test temp"
)

###
# [1st] How acc10 vs acc20 are responding at test temp 5, 10, 15, 20, 25, 32, 35
# Note: Ref is acclimation.temp=10 and test.temp=5
# At test temp = 5,
res_test5 <- results(dds_filt, name = "acclimation.temp_20_vs_10")
deg_test5 <- rownames(res_test5)[!is.na(res_test5$padj) & res_test5$padj < 0.05]
# At test temp = 10,
res_test10 <- results(dds_filt,
  contrast = list(c("acclimation.temp_20_vs_10","acclimation.temp20.test.temp10")))
deg_test10 <- rownames(res_test10)[!is.na(res_test10$padj) & res_test10$padj < 0.05]
# At test temp = 15,
res_test15 <- results(dds_filt,
  contrast = list(c("acclimation.temp_20_vs_10","acclimation.temp20.test.temp15")))
deg_test15 <- rownames(res_test15)[!is.na(res_test15$padj) & res_test15$padj < 0.05]
# At test temp = 20,
res_test20 <- results(dds_filt,
  contrast = list(c("acclimation.temp_20_vs_10","acclimation.temp20.test.temp20")))
deg_test20 <- rownames(res_test20)[!is.na(res_test20$padj) & res_test20$padj < 0.05]
# At test temp = 25,
res_test25 <- results(dds_filt,
  contrast = list(c("acclimation.temp_20_vs_10","acclimation.temp20.test.temp25")))
deg_test25 <- rownames(res_test25)[!is.na(res_test25$padj) & res_test25$padj < 0.05]
# At test temp = 30,
res_test30 <- results(dds_filt,
  contrast = list(c("acclimation.temp_20_vs_10","acclimation.temp20.test.temp30")))
deg_test30 <- rownames(res_test30)[!is.na(res_test30$padj) & res_test30$padj < 0.05]
# At test temp = 32,
res_test32 <- results(dds_filt,
  contrast = list(c("acclimation.temp_20_vs_10","acclimation.temp20.test.temp32")))
deg_test32 <- rownames(res_test32)[!is.na(res_test32$padj) & res_test32$padj < 0.05]
# At test temp = 35,
res_test35 <- results(dds_filt,
  contrast = list(c("acclimation.temp_20_vs_10","acclimation.temp20.test.temp35")))
deg_test35 <- rownames(res_test35)[!is.na(res_test35$padj) & res_test35$padj < 0.05]

# Venn diagram for 5, 20, 35
library(VennDiagram)
library(grid)
venn.plot <- venn.diagram(
  x = list(
    T5 = deg_test5,
    T20 = deg_test20,
    T35 = deg_test35
  ),
  filename = NULL,
  fill = c("skyblue","orange","tomato"),
  alpha = 0.5,
  cex = 1.4,
  cat.cex = 1.3
)
grid.newpage()
grid.draw(venn.plot)

# UpSet plot for all test temps
library(UpSetR)
gene_sets <- list(
  T5  = deg_test5,
  T10 = deg_test10,
  T15 = deg_test15,
  T20 = deg_test20,
  T25 = deg_test25,
  T30 = deg_test30,
  T32 = deg_test32,
  T35 = deg_test35
)
upset_data <- fromList(gene_sets)
upset(
  upset_data,
  order.by = "freq",
  sets.bar.color = "steelblue",
  main.bar.color = "black",
  keep.order = TRUE
)
# ComplexUpset plot for all test temps
library(ComplexUpset)
ComplexUpset::upset(
  upset_data,
  intersect = names(gene_sets)
)

# [2nd] How acc10 and acc20 are responding at change in temp -5, +0, +5, +10, +15
# 10at5 v 20at15 (-5)
res_10at5_20at15 <- results(
  dds_filt,
  contrast=list(c(
    "acclimation.temp_20_vs_10",
    "test.temp_15_vs_5",
    "acclimation.temp20.test.temp15"
  ))
)
deg_10at5_20at15 <- rownames(res_10at5_20at15)[!is.na(res_10at5_20at15$padj)
                                               & res_10at5_20at15$padj < 0.05]
# 10at10 v 20at20 (+0)
# First vector are positive terms, second vector are negative terms
# DESeq will automatically subtract second vector
res_10at10_20at20 <- results(
  dds_filt,
  contrast = list(
    c(
      "acclimation.temp_20_vs_10",
      "test.temp_20_vs_5",
      "acclimation.temp20.test.temp20"
    ),
    c(
      "test.temp_10_vs_5"
    )
  )
)
deg_10at10_20at20 <- rownames(res_10at10_20at20)[!is.na(res_10at10_20at20$padj)
                                                 & res_10at10_20at20$padj < 0.05]
# 10at15 v 20at25 (+5)
res_10at15_20at25 <- results(
  dds_filt,
  contrast = list(
    c(
      "acclimation.temp_20_vs_10",
      "test.temp_25_vs_5",
      "acclimation.temp20.test.temp25"
    ),
    c(
      "test.temp_15_vs_5"
    )
  )
)
deg_10at15_20at25 <- rownames(res_10at15_20at25)[!is.na(res_10at15_20at25$padj) 
                                                 & res_10at15_20at25$padj < 0.05]

# 10at20 v 20at30 (+10)
res_10at20_20at30 <- results(
  dds_filt,
  contrast = list(
    c(
      "acclimation.temp_20_vs_10",
      "test.temp_30_vs_5",
      "acclimation.temp20.test.temp30"
    ),
    c(
      "test.temp_20_vs_5"
    )
  )
)
deg_10at20_20at30 <- rownames(res_10at20_20at30)[!is.na(res_10at20_20at30$padj)
                                                 & res_10at20_20at30$padj < 0.05]
# 10at25 v 20at35 (+15)
res_10at25_20at35 <- results(
  dds_filt,
  contrast = list(
    c(
      "acclimation.temp_20_vs_10",
      "test.temp_35_vs_5",
      "acclimation.temp20.test.temp35"
    ),
    c(
      "test.temp_25_vs_5"
    )
  )
)
deg_10at25_20at35 <- rownames(res_10at25_20at35)[!is.na(res_10at25_20at35$padj)
                                                 & res_10at25_20at35$padj < 0.05]

deg_res_list <- list(
  diff_minus_5 = deg_10at5_20at15,
  diff_plus_0  = deg_10at10_20at20,
  diff_plus_5  = deg_10at15_20at25,
  diff_plus_10 = deg_10at20_20at30,
  diff_plus_15 = deg_10at25_20at35
)

# Venn diagram for -5, 0, +5, +10, +15
library(VennDiagram)
library(grid)
venn.plot <- venn.diagram(
  x = deg_res_list,
  filename = NULL,
  fill = c("skyblue","orange","lightgreen","tomato", "yellow"),
  alpha = 0.5,
  cex = 1.4,
  cat.cex = 1.3
)
grid.newpage()
grid.draw(venn.plot)

# Upset plot for -5, 0, +5, +10, +15
library(UpSetR)
upset_data <- fromList(deg_res_list)
UpSetR::upset(
  upset_data,
  order.by = "freq",
  sets.bar.color = "steelblue",
  main.bar.color = "black",
  keep.order = TRUE
)

# ComplexUpset plot for -5, 0, +5, +10, +15
library(ComplexUpset)
ComplexUpset::upset(
  upset_data,
  intersect = names(deg_res_list)
)

# [3rd] Effect of change in temp on acc10 and acc20 respectively
# Focus on acc10 ONLY - Relevel test temp to 10C so that ref is test.temp=10
dds_acc10 <- dds_filt
dds_acc10$test.temp <- relevel(
  dds_acc10$test.temp,
  ref = "10"
)
dds_acc10 <- DESeq(dds_acc10)

# Decrease temp by 5C for acc10
res_acc10_m5 <- results(
  dds_acc10,
  name = "test.temp_5_vs_10"
)
# Increase temp by 5C for acc10
res_acc10_p5 <- results(
  dds_acc10,
  name = "test.temp_15_vs_10"
)
# Increase temp by 10C for acc10
res_acc10_p10 <- results(
  dds_acc10,
  name = "test.temp_20_vs_10"
)
# Increase temp by 15C for acc10
res_acc10_p15 <- results(
  dds_acc10,
  name = "test.temp_25_vs_10"
)
# Increase temp by 20C for acc10
res_acc10_p20 <- results(
  dds_acc10,
  name = "test.temp_30_vs_10"
)
# Increase temp by 22C for acc10
res_acc10_p22 <- results(
  dds_acc10,
  name = "test.temp_32_vs_10"
)
# Increase temp by 25C for acc10
res_acc10_p25 <- results(
  dds_acc10,
  name = "test.temp_35_vs_10"
)
# Put into a list
acc10_response <- list(
  "-5"  = res_acc10_m5,
  "+5"  = res_acc10_p5,
  "+10" = res_acc10_p10,
  "+15" = res_acc10_p15,
  "+20" = res_acc10_p20,
  "+22" = res_acc10_p22,
  "+25" = res_acc10_p25
)
# Extract DEG counts
deg_counts_acc10 <- sapply(
  acc10_response,
  function(x){
    sum(
      !is.na(x$padj) &
        x$padj < 0.05 &
        abs(x$log2FoldChange) > 1
    )
  }
)
deg_counts_acc10

# Temp Response Curve for acc10

# Create function to count DEGs
count_DEGs <- function(res, padj_cutoff = 0.05, lfc_cutoff = 1){
  res <- as.data.frame(res)
  total <- sum(
    !is.na(res$padj) &
      res$padj < padj_cutoff &
      abs(res$log2FoldChange) > lfc_cutoff
  )
  up <- sum(
    !is.na(res$padj) &
      res$padj < padj_cutoff &
      res$log2FoldChange > lfc_cutoff
  )
  down <- sum(
    !is.na(res$padj) &
      res$padj < padj_cutoff &
      res$log2FoldChange < -lfc_cutoff
  )
  return(
    data.frame(
      Total = total,
      Up = up,
      Down = down
    )
  )
}
# Make results list
acc10_results <- list(
  "-5"  = res_acc10_m5,
  "+5"  = res_acc10_p5,
  "+10" = res_acc10_p10,
  "+15" = res_acc10_p15,
  "+20" = res_acc10_p20,
  "+22" = res_acc10_p22,
  "+25" = res_acc10_p25
)
# Extract DEG counts
deg_curve_acc10 <- do.call(
  rbind,
  lapply(names(acc10_results), function(x){
    out <- count_DEGs(acc10_results[[x]])
    data.frame(
      Delta_T = as.numeric(x),
      Total = out$Total,
      Up = out$Up,
      Down = out$Down
    )
  })
)
deg_curve_acc10

# Add 0C baseline point
deg_curve_acc10 <- rbind(
  data.frame(
    Delta_T = 0,
    Total = 0,
    Up = 0,
    Down = 0
  ),
  deg_curve_acc10
)
deg_curve_acc10 <- deg_curve_acc10[
  order(deg_curve_acc10$Delta_T),
]

# Convert to long format
library(tidyr)
library(dplyr)
deg_curve_long_acc10 <- deg_curve_acc10 %>%
  pivot_longer(
    cols = c(Total, Up, Down),
    names_to = "Response",
    values_to = "DEGs"
  )

# Plot response curve
library(ggplot2)
ggplot(
  deg_curve_long_acc10,
  aes(
    x = Delta_T,
    y = DEGs,
    group = Response,
    color = Response
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(0, 800)) +
  scale_color_manual(
    values = c(
      "Up" = "green",
      "Down" = "red",
      "Total" = "blue"
    )
  ) +
  scale_x_continuous(
    breaks = c(-5, 0, 5, 10, 15, 20, 22, 25),
    labels = c("-5", "0", "+5", "+10", "+15", "+20", "+22", "+25")
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Change in test temperature relative to acclimation (°C)",
    y = "Number of DEGs",
    title = "Transcriptomic response of 10°C-acclimated individuals"
  )

# Plot UpSet plot
get_DEGs <- function(res, padj_cutoff = 0.05, lfc_cutoff = 1){
  rownames(res)[
    !is.na(res$padj) &
      res$padj < padj_cutoff &
      abs(res$log2FoldChange) > lfc_cutoff
  ]
}
deg_sets_acc10 <- list(
  "-5"  = get_DEGs(res_acc10_m5),
  "+5"  = get_DEGs(res_acc10_p5),
  "+10" = get_DEGs(res_acc10_p10),
  "+15" = get_DEGs(res_acc10_p15),
  "+20" = get_DEGs(res_acc10_p20),
  "+22" = get_DEGs(res_acc10_p22),
  "+25" = get_DEGs(res_acc10_p25)
)
library(UpSetR)
upset_data <- fromList(deg_sets_acc10)
upset(
  upset_data,
  order.by = "freq",
  sets.bar.color = "steelblue",
  main.bar.color = "black",
  keep.order = TRUE,
  mainbar.y.label = "Number of shared DEGs",
  sets.x.label = "Number of DEGs per temperature change"
)
# ComplexUpset plot
library(ComplexUpset)
library(tidyverse)
# Get all unique DEG IDs
all_genes <- unique(unlist(deg_sets_acc10))
# Create presence/absence matrix
upset_df <- data.frame(
  gene = all_genes
)
for (set_name in names(deg_sets_acc10)) {
  upset_df[[set_name]] <- all_genes %in% deg_sets_acc10[[set_name]]
}
head(upset_df)

ComplexUpset::upset(
  upset_df,
  intersect = names(deg_sets_acc10),
  name = "Temperature change",
  sort_intersections_by = "cardinality",
  base_annotations = list(
    'Intersection size' = intersection_size()
  )
)

# Focus on acc20 ONLY - Relevel test temp to 20C so that ref is test.temp=20
dds_acc20 <- dds_filt
dds_acc20$test.temp <- relevel(
  dds_acc20$test.temp,
  ref="20"
)
dds_acc20 <- DESeq(dds_acc20)
# Decrease temp by 15C for acc20
res_acc20_m15 <- results(
  dds_acc20,
  name = "test.temp_5_vs_20"
)
# Decrease temp by 10C for acc20
res_acc20_m10 <- results(
  dds_acc20,
  name = "test.temp_10_vs_20"
)
# Decrease temp by 5C for acc20
res_acc20_m5 <- results(
  dds_acc20,
  name = "test.temp_15_vs_20"
)
# Increase temp by 5C for acc20
res_acc20_p5 <- results(
  dds_acc20,
  name = "test.temp_25_vs_20"
)
# Increase temp by 10C for acc20
res_acc20_p10 <- results(
  dds_acc20,
  name = "test.temp_30_vs_20"
)
# Increase temp by 12C for acc20
res_acc20_p12 <- results(
  dds_acc20,
  name = "test.temp_32_vs_20"
)
# Increase temp by 15C for acc20
res_acc20_p15 <- results(
  dds_acc20,
  name = "test.temp_35_vs_20"
)
# Put into a list
acc20_response <- list(
  "-15" = res_acc20_m15,
  "-10" = res_acc20_m10,
  "-5"  = res_acc20_m5,
  "+5"  = res_acc20_p5,
  "+10" = res_acc20_p10,
  "+12" = res_acc20_p12,
  "+15" = res_acc20_p15
)
# Extract DEG counts
deg_counts_acc20 <- sapply(
  acc20_response,
  function(x){
    sum(
      !is.na(x$padj) &
        x$padj < 0.05 &
        abs(x$log2FoldChange) > 1
    )
  }
)

# Temp Response Curve for acc20
# Create function to count DEGs
count_DEGs <- function(res, padj_cutoff = 0.05, lfc_cutoff = 1){
  res <- as.data.frame(res)
  total <- sum(
    !is.na(res$padj) &
      res$padj < padj_cutoff &
      abs(res$log2FoldChange) > lfc_cutoff
  )
  up <- sum(
    !is.na(res$padj) &
      res$padj < padj_cutoff &
      res$log2FoldChange > lfc_cutoff
  )
  down <- sum(
    !is.na(res$padj) &
      res$padj < padj_cutoff &
      res$log2FoldChange < -lfc_cutoff
  )
  return(
    data.frame(
      Total = total,
      Up = up,
      Down = down
    )
  )
}
# Make a results list
acc20_results <- list(
  "-15" = res_acc20_m15,
  "-10" = res_acc20_m10,
  "-5"  = res_acc20_m5,
  "+5"  = res_acc20_p5,
  "+10" = res_acc20_p10,
  "+12" = res_acc20_p12,
  "+15" = res_acc20_p15
)

# Extract DEG counts
deg_curve_acc20 <- do.call(
  rbind,
  lapply(names(acc20_results), function(x){
    
    out <- count_DEGs(acc20_results[[x]])
    
    data.frame(
      Delta_T = as.numeric(x),
      Total = out$Total,
      Up = out$Up,
      Down = out$Down
    )
  })
)

# Add 0C baseline point (Acc20 at test temperature 20C)
deg_curve_acc20 <- rbind(
  data.frame(
    Delta_T = 0,
    Total = 0,
    Up = 0,
    Down = 0
  ),
  deg_curve_acc20
)
# Sort by temperature change
deg_curve_acc20 <- deg_curve_acc20[
  order(deg_curve_acc20$Delta_T),
]
deg_curve_acc20

# Convert to long format for plotting
library(tidyr)
library(dplyr)
deg_curve_long <- deg_curve_acc20 %>%
  pivot_longer(
    cols = c(Total, Up, Down),
    names_to = "Response",
    values_to = "DEGs"
  )

# Plot the response curve
library(ggplot2)
ggplot(
  deg_curve_long,
  aes(
    x = Delta_T,
    y = DEGs,
    group = Response,
    color = Response
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(0, 800)) +
  scale_color_manual(
    values = c(
      "Up" = "green",
      "Down" = "red",
      "Total" = "blue"
    )
  ) +
  scale_x_continuous(
    breaks = c(-15, -10, -5, 0, 5, 10, 12, 15),
    labels = c("-15", "-10", "-5", "0", "+5", "+10", "+12", "+15")
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Change in test temperature relative to acclimation (°C)",
    y = "Number of DEGs",
    title = "Transcriptomic response of 20°C-acclimated individuals"
  )

# ComplexUpset plot
library(ComplexUpset)
library(tidyverse)
acc20_deg_lists <- lapply(acc20_results, function(res) {
  rownames(res)[which(res$padj < 0.05)]
})
all_genes <- unique(unlist(acc20_deg_lists))
acc20_upset_df <- data.frame(
  gene = all_genes
)
for (contrast in names(acc20_deg_lists)) {
  acc20_upset_df[[contrast]] <- all_genes %in% acc20_deg_lists[[contrast]]
}
head(acc20_upset_df)

ComplexUpset::upset(
  acc20_upset_df,
  names(acc20_deg_lists),
  min_size = 2,
  width_ratio = 0.2,
  base_annotations = list(
    'Intersection size' = intersection_size()
  )
)

#####
# Make a table of DEG count for each pairwise comparison
# 21 Oct 2026: added min absolute log2foldchange of 1
get_DEGs <- function(res, padj_cutoff = 0.05, lfc_cutoff = 1) {
  res <- as.data.frame(res)
  rownames(res)[
    !is.na(res$padj) &
      res$padj < padj_cutoff &
      !is.na(res$log2FoldChange) &
      abs(res$log2FoldChange) >= lfc_cutoff
  ]
}
deg_sets_acc10 <- lapply(
  acc10_results,
  get_DEGs
)
deg_sets_acc20 <- lapply(
  acc20_results,
  get_DEGs
)
sapply(deg_sets_acc10, length)
sapply(deg_sets_acc20, length)

# Make a dot plot of shared_matrix
library(ggplot2)
shared_df <- as.data.frame(as.table(shared_matrix))
colnames(shared_df) <- c("acc10", "acc20", "shared_genes")
shared_df$shared_genes <- as.numeric(shared_df$shared_genes)
# Add original DEG count for each acc comparison
acc10_deg_counts <- sapply(deg_sets_acc10, length)
shared_df$original_acc10_DEGs <- acc10_deg_counts[
  as.character(shared_df$acc10)]
acc20_deg_counts <- sapply(deg_sets_acc20, length)
shared_df$original_acc20_DEGs <- acc20_deg_counts[
  as.character(shared_df$acc20)
]
# Set order and labels for y-axis (acc10)
shared_df$acc10 <- factor(
  shared_df$acc10,
  levels = rownames(shared_matrix),
  labels = c(
    "-15\n(5)",
    "-10\n(10)",
    "-5\n(15)",
    "+5\n(25)",
    "+10\n(30)",
    "+12\n(32)",
    "+15\n(35)"
  )
)
# Set order and labels for x-axis (acc20)
shared_df$acc20 <- factor(
  shared_df$acc20,
  levels = colnames(shared_matrix),
  labels = c(
    "-15\n(5)",
    "-10\n(10)",
    "-5\n(15)",
    "+5\n(25)",
    "+10\n(30)",
    "+12\n(32)",
    "+15\n(35)"
  )
)
# Plot
ggplot(
  shared_df,
  aes(
    x = acc20,
    y = acc10,
    size = shared_genes,
    color = shared_genes
  )
) +
  geom_point(alpha = 0.8) +
  geom_text(
    aes(
      label = paste0(
        original_acc10_DEGs,
        "/",
        original_acc20_DEGs
      )
    ),
    hjust = -0.5,
    size = 2,
    color = "black"
  ) +
  scale_size_continuous(
    name = "Shared DEGs"
  ) +
  scale_color_gradient(
    name = "Shared DEGs",
    low = "lightgrey",
    high = "darkblue"
  ) +
  labs(
    x = "20°C-acclimated: ΔT (test temperature, °C)",
    y = "10°C-acclimated: ΔT (test temperature, °C)"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    legend.title = element_text(size = 11)
  ) +
  coord_cartesian(clip = "off")

