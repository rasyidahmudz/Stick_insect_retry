### WGCNA SCRIPT ###
setwd("~/Desktop/Stick_insect")
library(WGCNA)
allowWGCNAThreads()

# Load dataset [quantified and normalized]
dds_filt <- readRDS("~/Desktop/Stick_insect/rds/dds_filt.rds")
sampledata <- readRDS("~/Desktop/Stick_insect/rds/sampledata.rds")

# Extract expression matrix
library(SummarizedExperiment)
library(DESeq2)
vst_mat <- assay(vst(dds_filt))
# Samples should be rows for WGCNA
datExpr0 <- t(vst_mat)

# Trait data
datTraits <- sampledata[, c("acclimation.temp", "test.temp", "mass")]
# Save sample IDs from sampledata row names
sample_ids <- rownames(datTraits)
# Make sure categorical variables are factors
datTraits$acclimation.temp <- factor(datTraits$acclimation.temp)
datTraits$test.temp <- factor(datTraits$test.temp)
datTraits$mass <- as.numeric(datTraits$mass)

# Convert to numeric design matrix (safe for WGCNA) - TRIAL!!!
datTraits <- model.matrix(
  ~ -1 + acclimation.temp:test.temp + mass,
  data = datTraits
)
# Restore sample IDs as rownames
rownames(datTraits) <- sample_ids

# Reorder rows
datTraits <- datTraits[rownames(datExpr0), , drop = FALSE]
all(rownames(datExpr0) == rownames(datTraits)) # sanity check!

# Pick soft threshold
sft <- pickSoftThreshold(datExpr0,
                         powerVector = 1:20,
                         networkType = "signed",
                         verbose = 5)
# Plot
par(mfrow = c(1, 2))
# Left: Scale-free topology fit
plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Fit Index",
  type = "b",
  main = "Scale independence"
)
# Right: Mean connectivity
plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "b",
  main = "Mean connectivity"
)
dev.off()

# Network construction
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(datExpr0,                # <= input here
                          
                          # == Adjacency Function ==
                          power = 9,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

netwk$colors[netwk$blockGenes[[1]]]
table(netwk$colors)

# Relate module assignments to treatment groups
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df[1:5,]

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(datExpr0, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)

# Module-trait correlations
moduleTraitCor <- cor(MEs0, datTraits, use = "p")

moduleTraitP <- corPvalueStudent(
  moduleTraitCor,
  nSamples = nrow(datExpr0)
)

x_labels <- colnames(moduleTraitCor)
y_labels <- gsub("^ME", "", rownames(moduleTraitCor))

# Desired order of test temperatures
test_temps <- c(5, 10, 15, 20, 25, 30, 32, 35)

# Create the desired order of x-axis labels
desired_order <- c(
  "mass",
  paste0("acclimation.temp10:test.temp", test_temps),
  paste0("acclimation.temp20:test.temp", test_temps)
)

# Reorder the heatmap correlation matrix
moduleTraitCor <- moduleTraitCor[, desired_order, drop = FALSE]

# Reorder the p-value matrix too, if used
moduleTraitP <- moduleTraitP[, desired_order, drop = FALSE]

# Update x-axis labels
x_labels <- colnames(moduleTraitCor)

# Plot heatmap
png(
  "newmap_oct29.png",
  width = 5000,
  height = 3000,
  res = 300
)

# Increase the left margin to make room for y-axis labels
par(mar = c(8, 12, 4, 2) + 0.1)

newmap <- labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = x_labels,
  yLabels = y_labels,
  ySymbols = y_labels,
  colors = blueWhiteRed(50),
  zlim = c(-0.6, 0.6),
  textMatrix = signif(moduleTraitCor, 2),
  setStdMargins = FALSE,
  cex.text = 0.7,
  cex.lab.x = 0.8,
  cex.lab.y = 0.8,
  xLabelsAngle = 25,
  main = "Module–Trait Relationships"
)
dev.off()


### Reduced dimensionality heatmap - CORRELATION ###

# Test temperatures in the desired order
test_temps <- c(5, 10, 15, 20, 25, 30, 32, 35)

# Create reduced matrix: 20C-acclimated minus 10C-acclimated
reducedCor <- sapply(test_temps, function(temp) {
  moduleTraitCor[, paste0("acclimation.temp20:test.temp", temp)] -
    moduleTraitCor[, paste0("acclimation.temp10:test.temp", temp)]
})

# Set column names
colnames(reducedCor) <- paste0("Test temp ", test_temps, "°C")

# Keep module names as row names
rownames(reducedCor) <- rownames(moduleTraitCor)

png(
  "reduced_newmap_oct29.png",
  width = 5000,
  height = 3000,
  res = 300
)

par(mar = c(7, 12, 4, 2) + 0.1)

reducedMap <- labeledHeatmap(
  Matrix = reducedCor,
  xLabels = colnames(reducedCor),
  yLabels = gsub("^ME", "", rownames(reducedCor)),
  ySymbols = gsub("^ME", "", rownames(reducedCor)),
  colors = blueWhiteRed(50),
  zlim = c(-0.6, 0.6),
  textMatrix = signif(reducedCor, 2),
  setStdMargins = FALSE,
  cex.text = 0.8,
  cex.lab.x = 0.9,
  cex.lab.y = 0.8,
  main = "Difference in Module-Trait Correlations\n(20°C minus 10°C acclimation)"
)

dev.off()

### COUNT NUMBER OF GENES ###
# Total genes included in the network
length(netwk$colors)

# Total genes assigned to non-grey modules
sum(labels2colors(netwk$colors) != "grey")

# Create a sorted summary table
module_counts <- as.data.frame(
  table(Module = labels2colors(netwk$colors))
)

colnames(module_counts) <- c("Module", "Number_of_genes")

module_counts <- module_counts[
  order(-module_counts$Number_of_genes),
]

rownames(module_counts) <- NULL
module_counts

# Plot the no. of genes in modules respectively, excluding grey module
module_counts_plot <- subset(module_counts, Module != "grey")
library(ggplot2)
# Order modules by gene count
module_counts$Module <- reorder(
  module_counts$Module,
  module_counts$Number_of_genes
)

ggplot(module_counts, aes(x = Module, y = Number_of_genes)) +
  geom_col() +
  coord_flip() +
  geom_text(
    aes(label = Number_of_genes),
    hjust = -0.15,
    size = 3
  ) +
  labs(
    title = "Number of Genes per WGCNA Module",
    x = "WGCNA Module",
    y = "Number of Genes"
  ) +
  expand_limits(y = max(module_counts$Number_of_genes) * 1.12) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 9)
  )


### RETAINED MODULES: FILTER, NUMBER, SUMMARISE AND PLOT ###

library(dplyr)
library(ggplot2)
library(readr)
library(rtracklayer)
library(gridExtra)
library(grid)

## 1. Filter modules using the existing correlation matrix ##

# Select acclimation x test-temperature columns
temp_cols <- grep(
  "acclimation.temp",
  colnames(moduleTraitCor),
  value = TRUE
)

# Calculate maximum absolute correlation for each module
max_abs_cor <- apply(
  abs(moduleTraitCor[, temp_cols, drop = FALSE]),
  1,
  max,
  na.rm = TRUE
)

# Retain modules with at least one |r| >= 0.3, excluding grey
modules_to_keep <- names(max_abs_cor)[
  max_abs_cor >= 0.3 &
    sub("^ME", "", names(max_abs_cor)) != "grey"
]

# Confirm the expected number of retained modules
stopifnot(length(modules_to_keep) == 18)


## 2. Get gene assignments and count genes per retained module ##

# Get module colour assignments from the current WGCNA result
module_colors <- labels2colors(netwk$colors)

# Use the gene IDs associated with the WGCNA assignments
gene_ids <- names(netwk$colors)

if (is.null(gene_ids)) {
  gene_ids <- colnames(datExpr0)
}

gene_modules <- data.frame(
  gene_id = gene_ids,
  Module_Color = as.character(module_colors)
)

# Count genes in retained modules only
module_counts <- gene_modules %>%
  filter(Module_Color %in% sub("^ME", "", modules_to_keep)) %>%
  count(Module_Color, name = "Genes_in_Module") %>%
  arrange(desc(Genes_in_Module), Module_Color)

# Assign module numbers in descending order of gene count
retained_module_key <- module_counts %>%
  mutate(
    Module_Number = sprintf("%02d", row_number()),
    Module_ID = paste0("Module_", Module_Number)
  ) %>%
  select(Module_ID, Module_Number, Module_Color, Genes_in_Module)

stopifnot(nrow(retained_module_key) == 18)

# Total number of genes in retained modules
total_retained_genes <- sum(retained_module_key$Genes_in_Module)

retained_module_key
total_retained_genes


## 3. Create gene-to-module assignments ##

# Add module numbers to retained genes.
# Excluded and grey modules remain unnumbered.
gene_modules <- gene_modules %>%
  left_join(
    retained_module_key,
    by = "Module_Color"
  ) %>%
  mutate(
    Module_Type = case_when(
      Module_Color == "grey" ~ "unassigned",
      !is.na(Module_ID) ~ "retained",
      TRUE ~ "excluded"
    )
  )


## 4. Export simplified module summary ##

Module_Summary_Simplified <- retained_module_key %>%
  select(
    Module_ID,
    Module_Number,
    Module_Color,
    Genes_in_Module
  )

write_csv(
  Module_Summary_Simplified,
  "clitarchus_Module_Summary_Simplified.csv"
)


## 5. Export annotated module summary ##

# Import GTF annotation
gtf <- rtracklayer::import("cli_annevo_gemoma_combined.gtf")
gtf_df <- as.data.frame(gtf)

# Keep gene-level records where available
if ("type" %in% colnames(gtf_df) &&
    any(gtf_df$type == "gene", na.rm = TRUE)) {
  gtf_df <- gtf_df %>% filter(type == "gene")
}

# Select available annotation fields
annotation_fields <- intersect(
  c("gene_id", "gene_name", "product", "description"),
  colnames(gtf_df)
)

gene_annotation <- gtf_df %>%
  select(all_of(annotation_fields)) %>%
  filter(!is.na(gene_id)) %>%
  distinct(gene_id, .keep_all = TRUE)

# Join annotations to module assignments
Module_Summary <- gene_modules %>%
  left_join(gene_annotation, by = "gene_id") %>%
  arrange(
    factor(Module_Type, levels = c("retained", "excluded", "unassigned")),
    Module_Number,
    gene_id
  )

write_csv(
  Module_Summary,
  "clitarchus_Module_Summary_Annotated.csv"
)


## 6. Plot gene counts for the 18 retained modules ##

module_counts_plot <- retained_module_key %>%
  mutate(
    Module_ID = factor(
      Module_ID,
      levels = rev(Module_ID)
    )
  )

gene_barplot <- ggplot(
  module_counts_plot,
  aes(x = Module_ID, y = Genes_in_Module)
) +
  geom_col() +
  geom_text(
    aes(label = Genes_in_Module),
    hjust = -0.15,
    size = 3
  ) +
  coord_flip() +
  labs(
    title = "Number of Genes in Retained WGCNA Modules",
    subtitle = paste(
      "18 modules |",
      format(total_retained_genes, big.mark = ","),
      "genes"
    ),
    x = "WGCNA Module",
    y = "Number of Genes"
  ) +
  expand_limits(y = max(module_counts_plot$Genes_in_Module) * 1.12) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(
  "retained_WGCNA_module_gene_counts.png",
  plot = gene_barplot,
  width = 7,
  height = 6,
  dpi = 300,
  bg = "white"
)


## 7. Export module summary table as a figure ##

module_table <- retained_module_key %>%
  select(
    Module_ID,
    Module_Color,
    Genes_in_Module
  )

table_plot <- tableGrob(
  module_table,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 10
  )
)

png(
  "retained_WGCNA_module_summary.png",
  width = 1800,
  height = 1200,
  res = 200
)

grid.draw(table_plot)
dev.off()


## 8. Create numbered heatmap using the existing filtered correlation matrix ##
## Rebuild heatmap using all 18 retained modules ##

# Module order is already defined by retained_module_key
heatmap_order <- paste0(
  "ME",
  retained_module_key$Module_Color
)

# Check that all 18 modules exist in the original correlation matrix
setdiff(heatmap_order, rownames(moduleTraitCor))

# Extract all 18 modules in the correct order
filteredCor <- moduleTraitCor[
  heatmap_order,
  ,
  drop = FALSE
]

# Create matching numbered labels
filtered_y_labels <- paste0(
  retained_module_key$Module_ID,
  " (",
  retained_module_key$Genes_in_Module,
  ")"
)

# Confirm
dim(filteredCor)
length(filtered_y_labels)

# Plot numbered heatmap
png(
  "retained_WGCNA_module_trait_heatmap.png",
  width = 5000,
  height = 3000,
  res = 300
)

par(mar = c(8, 12, 4, 2) + 0.1)

labeledHeatmap(
  Matrix = filteredCor,
  xLabels = x_labels,
  yLabels = filtered_y_labels,
  ySymbols = filtered_y_labels,
  colors = blueWhiteRed(50),
  zlim = c(-0.5, 0.5),
  textMatrix = signif(filteredCor, 2),
  setStdMargins = FALSE,
  cex.text = 0.7,
  cex.lab.x = 0.8,
  cex.lab.y = 0.8,
  xLabelsAngle = 25,
  main = "Module-Trait Relationships: 18 Retained Modules"
)

dev.off()

### SAVE OBJECTS ###
save(
  netwk,
  datExpr0,
  datTraits,
  moduleTraitCor,
  moduleTraitP,
  filteredCor,
  retained_module_key,
  gene_modules,
  file = "~/Desktop/Stick_insect/retry_obj/WGCNA_objects.RData"
)

### RESTORE OBJECTS ###
load("~/Desktop/Stick_insect/retry_obj/WGCNA_objects.RData")


