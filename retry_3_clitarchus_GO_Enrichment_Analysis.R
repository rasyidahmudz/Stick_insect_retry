################################################################################
#
# Title: WGCNA Module GO Enrichment Analysis
#        Colorado Potato Beetle (CPB) Diapause Gene Co-expression Networks
#
# Author: Alex Torson
# Institution: USDA-ARS
# Email: Alex.Torson@usda.gov
#
# Description: Gene Ontology enrichment analysis on WGCNA modules across
#              Biological Process, Molecular Function, and Cellular Component
#              ontologies to identify functional themes within co-expression
#              modules.
#
# Dependencies: Must run 09_core_wgcna_analysis.R first
#
# Input files:
#   - 09_core_results/Module_Summary.csv
#   - Annotation/cpb_diapause_go_omicsbox_export_one_seq_per_row.txt
#
# Output files:
#   - 12_go_enrichment/module_enrichments/*.csv (individual module results)
#   - 12_go_enrichment/summary_tables/*.csv (combined results)
#   - 12_go_enrichment/plots/*.png (visualization plots)
#
################################################################################

# Package Loading and Setup -----------------------------------------------

library(tidyverse)
library(topGO)
library(rtracklayer)

# Working Directory and Output Setup --------------------------------------

setwd("~/Desktop/Stick_insect/Stick_insect_retry/enrichment_analyses")

go_dirs <- c(
  "go/module_enrichments",
  "go/summary_tables",
  "go/plots"
)

walk(go_dirs, ~ if (!dir.exists(.x)) dir.create(.x, recursive = TRUE))

# Standardized Color Palette and Theme ------------------------------------

# Diapause status colors (primary scheme)
diapause_colors <- c("Diapause" = "#0072B2", "NonDiapause" = "#E69F00")

# Sex colors (assigned from variable clusters)
sex_colors <- c("Male" = "#00AFBB", "Female" = "#009E73")

# K-means cluster colors (for PCA variable plots)
cluster_colors <- c("#E69F00", "#00AFBB", "#009E73", "#0072B2")

# Significance colors for MA plots
sig_colors <- c("TRUE" = "#D55E00", "FALSE" = "black")

# Standardized theme
theme_cpb <- theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'plain', size = 12),
    axis.text = element_text(face = "plain", size = 10),
    axis.title = element_text(face = "plain", size = 10),
    axis.title.x = element_text(margin = margin(t = 8, r = 20, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 6, b = 0, l = 0)),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# Data Loading and Preparation --------------------------------------------

# Load WGCNA module assignments
module_assignments <- read_csv("~/Desktop/Stick_insect/clitarchus_Module_Summary_Annotated.csv")

# Load GTF annotations
cli_annot <- rtracklayer::import(
  "~/Desktop/Stick_insect/cli_annevo_gemoma_combined.gtf"
) %>%
  as.data.frame()

# Create transcript-to-gene mapping
transcript_to_gene <- cli_annot %>%
  dplyr::select(transcript_id, gene_id) %>%
  filter(!is.na(transcript_id), !is.na(gene_id)) %>%
  distinct()

# Load BLAST annotations
blast_annot <- read.csv(
  "~/Desktop/Stick_insect/clitarchus_annotation.csv"
)

# Create gene annotation lookup
gene_annotation_lookup <- blast_annot %>%
  dplyr::select(
    gene_id = qry_gene_id,
    gene_annotation
  ) %>%
  distinct(gene_id, .keep_all = TRUE)

# Load original GO annotations file
go_raw <- read.delim(
  "~/Desktop/Stick_insect/clitarchus_seq_GO.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Convert comma-separated GO IDs into individual rows
go_annotations <- go_raw %>%
  filter(
    !is.na(Annotation.GO.ID),
    Annotation.GO.ID != ""
  ) %>%
  separate_rows(Annotation.GO.ID, sep = ",") %>%
  mutate(
    go_id = trimws(Annotation.GO.ID)
  ) %>%
  dplyr::select(
    Sequence.Name,
    go_id
  ) %>%
  distinct()

# Map transcript IDs to gene IDs using your GTF mapping
transcript_to_gene <- cli_annot %>%
  dplyr::select(transcript_id, gene_id) %>%
  distinct()

go_annotations <- go_annotations %>%
  inner_join(
    transcript_to_gene,
    by = c("Sequence.Name" = "transcript_id")
  ) %>%
  dplyr::select(gene_id, go_id) %>%
  distinct()

# Add gene annotations to GO annotation table
go_annotations_annotated <- go_annotations %>%
  left_join(
    gene_annotation_lookup,
    by = "gene_id"
  )

# Check
head(go_annotations_annotated)

# -------------------------------------------------------
# DEFINE EXPRESSED GENE UNIVERSE
# -------------------------------------------------------
expressed_gene_ids <- module_assignments %>%
  pull(gene_id) %>%
  unique()

# Confirm expected number of expressed genes
length(expressed_gene_ids)
stopifnot(length(expressed_gene_ids) == 19252)

# Retain GO annotations for expressed genes only
go_annotations <- go_annotations %>%
  filter(gene_id %in% expressed_gene_ids) %>%
  distinct(gene_id, go_id)

# Create gene-to-GO mapping list
gene2GO <- split(
  go_annotations$go_id,
  go_annotations$gene_id
)

# Remove duplicate GO terms per gene
gene2GO <- lapply(gene2GO, unique)

# Define universe
geneUniverse <- expressed_gene_ids

# Checks!
cat("Total expressed genes:",
    length(geneUniverse), "\n")
cat("Expressed genes with GO annotations:",
    n_distinct(go_annotations$gene_id), "\n")
cat("Unique GO terms:",
    n_distinct(go_annotations$go_id), "\n")

# -------------------------------------------------------
# DEFINE MODULES FOR GO ENRICHMENT
# -------------------------------------------------------

# Filter to genes with GO annotations
module_assignments_go <- module_assignments %>%
  filter(gene_id %in% go_annotations$gene_id)

# Select retained modules only
modules_to_analyze <- module_assignments_go %>%
  filter(Module_Type == "retained") %>%
  group_by(Module_ID) %>%
  summarise(
    genes = list(gene_id),
    n_genes = n_distinct(gene_id),
    .groups = "drop"
  ) %>%
  filter(n_genes >= 10)

# Check modules
modules_to_analyze |> 
  dplyr::select(Module_ID, n_genes)


# -------------------------------------------------------
# GO ENRICHMENT FUNCTION
# -------------------------------------------------------

perform_module_go_enrichment <- function(
    module_id,
    module_genes,
    ontology = "BP") {
  
  # Create binary gene list using the full expressed gene universe
  geneList <- factor(
    as.integer(geneUniverse %in% module_genes)
  )
  
  names(geneList) <- geneUniverse
  
  # Check number of genes in module
  genes_in_module <- sum(geneList == 1)
  
  if (genes_in_module < 5) {
    return(NULL)
  }
  
  # Create topGO data object
  GOdata <- new(
    "topGOdata",
    description = paste(module_id, ontology),
    ontology = ontology,
    allGenes = geneList,
    annot = annFUN.gene2GO,
    gene2GO = gene2GO,
    nodeSize = 5
  )
  
  # Run parent-child Fisher's exact test
  resultFisher <- runTest(
    GOdata,
    algorithm = "parentchild",
    statistic = "fisher"
  )
  
  # Get all tested GO terms
  allGO <- usedGO(GOdata)
  
  if (length(allGO) == 0) {
    return(NULL)
  }
  
  # Extract numerical p-values directly
  go_scores <- topGO::score(resultFisher)
  
  # Generate GO term descriptions and statistics
  allRes <- GenTable(
    GOdata,
    p_value = resultFisher,
    topNodes = length(allGO),
    numChar = 1000
  )
  
  # Replace formatted p-values with numerical values
  allRes <- allRes %>%
    mutate(
      p_value = as.numeric(go_scores[GO.ID])
    ) %>%
    filter(!is.na(p_value))
  
  # Calculate BH-adjusted p-values
  allRes <- allRes %>%
    mutate(
      p.adjust = p.adjust(p_value, method = "BH"),
      Significant = as.numeric(Significant),
      Expected = as.numeric(Expected),
      module_id = module_id,
      ontology = ontology,
      enrichment_ratio = ifelse(
        Expected > 0,
        Significant / Expected,
        NA_real_
      )
    ) %>%
    arrange(p.adjust)
  
  # Extract genes contributing to each GO term
  allRes$genes <- vapply(
    allRes$GO.ID,
    function(go_term) {
      
      genes_in_term <- genesInTerm(
        GOdata,
        go_term
      )[[1]]
      
      module_genes_in_term <- intersect(
        genes_in_term,
        module_genes
      )
      
      paste(module_genes_in_term, collapse = "; ")
      
    },
    character(1)
  )
  
  # Extract functional annotations for contributing genes
  allRes$gene_annotations <- vapply(
    allRes$genes,
    function(gene_list) {
      
      genes <- strsplit(gene_list, "; ", fixed = TRUE)[[1]]
      genes <- genes[nzchar(genes)]
      
      annotations <- gene_annotation_lookup$gene_annotation[
        match(genes, gene_annotation_lookup$gene_id)
      ]
      
      annotations[is.na(annotations) | annotations == ""] <- "No annotation"
      
      paste(annotations, collapse = "; ")
      
    },
    character(1)
  )
  
  # Combine gene IDs and functional annotations
  allRes$gene_details <- vapply(
    allRes$genes,
    function(gene_list) {
      
      genes <- strsplit(gene_list, "; ", fixed = TRUE)[[1]]
      genes <- genes[nzchar(genes)]
      
      annotations <- gene_annotation_lookup$gene_annotation[
        match(genes, gene_annotation_lookup$gene_id)
      ]
      
      annotations[is.na(annotations) | annotations == ""] <- "No annotation"
      
      details <- paste0(
        genes,
        " (",
        annotations,
        ")"
      )
      
      paste(details, collapse = "; ")
      
    },
    character(1)
  )
  
  # Retain FDR-significant GO terms
  allRes <- allRes %>%
    filter(p.adjust <= 0.05)
  
  return(allRes)
}

# Run GO Enrichment -------------------------------------------------------

all_enrichment_results <- list()
go_errors <- list()
ontologies <- c("BP", "MF", "CC")

for (ont in ontologies) {
  
  ont_results_list <- list()
  
  for (i in seq_len(nrow(modules_to_analyze))) {
    
    module_id <- modules_to_analyze$Module_ID[i]
    module_genes <- modules_to_analyze$genes[[i]]
    
    message("Running ", module_id, " - ", ont)
    
    module_result <- tryCatch({
      
      perform_module_go_enrichment(
        module_id,
        module_genes,
        ont
      )
      
    }, error = function(e) {
      
      go_errors[[paste(module_id, ont, sep = "_")]] <<- conditionMessage(e)
      
      message(
        "ERROR: ", module_id, " - ", ont,
        ": ", conditionMessage(e)
      )
      
      NULL
    })
    
    if (!is.null(module_result) && nrow(module_result) > 0) {
      ont_results_list[[module_id]] <- module_result
    }
  }
  
  if (length(ont_results_list) > 0) {
    all_enrichment_results[[ont]] <- bind_rows(
      ont_results_list
    )
  }
}

# Save Individual Module Results ------------------------------------------

for (i in 1:nrow(modules_to_analyze)) {
  module_id <- modules_to_analyze$Module_ID[i]
  module_genes <- modules_to_analyze$genes[[i]]
  
  # Create module directory
  module_dir <- file.path("./go/module_enrichments", module_id)
  if (!dir.exists(module_dir)) dir.create(module_dir, recursive = TRUE)
  
  # Save results for each ontology
  for (ont in ontologies) {
    if (ont %in% names(all_enrichment_results)) {
      existing_result <- all_enrichment_results[[ont]] %>%
        filter(module_id == !!module_id)
      
      if (nrow(existing_result) > 0) {
        write_csv(existing_result, file.path(module_dir, paste0(ont, "_enrichment.csv")))
      }
    }
  }
}

# -------------------------------------------------------
# CREATE OVERALL MODULE ENRICHMENT SUMMARIES
# -------------------------------------------------------

# Create all module-ontology combinations
module_enrichment_summary <- expand_grid(
  module_id = modules_to_analyze$Module_ID,
  ontology_type = ontologies
)

# Prepare enrichment counts
if (length(all_enrichment_results) > 0) {
  
  enrichment_counts <- bind_rows(
    all_enrichment_results
  ) %>%
    group_by(module_id, ontology) %>%
    summarise(
      n_significant_terms = n(),
      top_p_value = min(p_value, na.rm = TRUE),
      top_p_adjust = min(p.adjust, na.rm = TRUE),
      avg_enrichment_ratio = mean(
        enrichment_ratio,
        na.rm = TRUE
      ),
      top_GO_term = Term[which.min(p.adjust)],
      .groups = "drop"
    ) %>%
    rename(ontology_type = ontology)
  
} else {
  
  # Create empty summary with required columns
  enrichment_counts <- tibble(
    module_id = character(),
    ontology_type = character(),
    n_significant_terms = integer(),
    top_p_value = numeric(),
    top_p_adjust = numeric(),
    avg_enrichment_ratio = numeric(),
    top_GO_term = character()
  )
  
}

# Join enrichment counts to all modules
module_enrichment_summary <- module_enrichment_summary %>%
  left_join(
    enrichment_counts,
    by = c("module_id", "ontology_type")
  ) %>%
  mutate(
    n_significant_terms = replace_na(
      n_significant_terms, 0L
    ),
    has_enrichment = n_significant_terms > 0
  )

# Save summary
write_csv(
  module_enrichment_summary,
  "./go/summary_tables/Module_GO_Enrichment_Summary.csv"
)

# Save combined enrichment results

if (length(all_enrichment_results) > 0) {
  
  combined_enrichments <- bind_rows(
    all_enrichment_results,
    .id = "ontology"
  )
  
  write_csv(
    combined_enrichments,
    "./go/summary_tables/All_Modules_GO_Enrichment.csv"
  )
  
} else {
  
  combined_enrichments <- tibble()
  
  message("No significant GO enrichment results found.")
  
}

# Module Enrichment Status

module_go_status <- module_enrichment_summary %>%
  group_by(module_id) %>%
  summarise(
    BP = any(ontology_type == "BP" & has_enrichment),
    MF = any(ontology_type == "MF" & has_enrichment),
    CC = any(ontology_type == "CC" & has_enrichment),
    any_GO_enrichment = any(has_enrichment),
    total_significant_terms = sum(n_significant_terms),
    .groups = "drop"
  )

write_csv(
  module_go_status,
  "./go/summary_tables/Module_GO_Enrichment_Status.csv"
)

# Check results!
print(module_enrichment_summary)
print(module_go_status)


# Create Visualization Plots ----------------------------------------------

if (length(all_enrichment_results) > 0 && nrow(module_enrichment_summary) > 0) {
  
  # Heatmap
  enrichment_heatmap_data <- module_enrichment_summary %>%
    filter(n_significant_terms > 0) %>%
    mutate(
      ontology_name = case_when(
        ontology_type == "BP" ~ "Biological Process",
        ontology_type == "MF" ~ "Molecular Function", 
        ontology_type == "CC" ~ "Cellular Component"
      )
    )
  
  if (nrow(enrichment_heatmap_data) > 0) {
    p_heatmap <- ggplot(enrichment_heatmap_data, 
                        aes(x = ontology_name, y = module_id, fill = n_significant_terms)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = n_significant_terms), color = "black", size = 3) +
      scale_fill_gradient(low = "white", high = sig_colors[["TRUE"]], name = "Significant\nTerms") +
      labs(
        x = "GO Ontology",
        y = "WGCNA Module"
      ) +
      theme_cpb +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
      )
    
    ggsave(
      file.path("./go/plots", "Module_GO_Enrichment_Heatmap.pdf"),
      p_heatmap, width = 8, height = 10, dpi = 300
    )
    
    ggsave(
      file.path("./go/plots", "Module_GO_Enrichment_Heatmap.png"),
      p_heatmap, width = 8, height = 10, dpi = 300
    )
  }
  
  # Module enrichment summary bar plot
  enrichment_summary_plot_data <- module_enrichment_summary %>%
    group_by(ontology_type) %>%
    summarise(
      total_modules = n_distinct(module_id),
      modules_with_enrichment = sum(n_significant_terms > 0),
      total_terms = sum(n_significant_terms),
      .groups = 'drop'
    ) %>%
    mutate(
      ontology_name = case_when(
        ontology_type == "BP" ~ "Biological Process",
        ontology_type == "MF" ~ "Molecular Function",
        ontology_type == "CC" ~ "Cellular Component"
      )
    )
  
  if (nrow(enrichment_summary_plot_data) > 0) {
    p_summary <- ggplot(enrichment_summary_plot_data, 
                        aes(x = ontology_name, y = total_terms, fill = ontology_type)) +
      geom_col(color = "black", linewidth = 0.5) +
      geom_text(aes(label = total_terms), vjust = -0.5, size = 4) +
      scale_fill_manual(values = cluster_colors[1:3]) +
      labs(
        x = "GO Ontology",
        y = "Total Significant Terms"
      ) +
      theme_cpb +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    ggsave(
      file.path("./go/plots", "GO_Enrichment_Summary_Barplot.pdf"),
      p_summary, width = 8, height = 6, dpi = 300
    )
    
    ggsave(
      file.path("./go/plots", "GO_Enrichment_Summary_Barplot.png"),
      p_summary, width = 8, height = 6, dpi = 300
    )
  }
}

# Analysis Summary ---------------------------------------------------------

if (length(all_enrichment_results) > 0) {
  for (ont in names(all_enrichment_results)) {
    ont_data <- all_enrichment_results[[ont]]
    cat(sprintf("%s: %d significant terms across %d modules\n", 
                ont, nrow(ont_data), length(unique(ont_data$module_id))))
  }
} else {
  cat("No significant enrichments found.\n")
}

# -------------------------------------------------------
# SAVE GO ANALYSIS SUMMARY
# -------------------------------------------------------

analysis_summary <- list(
  
  modules_analyzed = nrow(modules_to_analyze),
  
  modules_with_enrichments = if (
    length(all_enrichment_results) > 0
  ) {
    length(unique(
      bind_rows(all_enrichment_results)$module_id
    ))
  } else {
    0
  },
  
  total_expressed_genes = length(geneUniverse),
  
  genes_with_GO_annotations = n_distinct(
    go_annotations$gene_id
  ),
  
  total_GO_terms = n_distinct(
    go_annotations$go_id
  ),
  
  enrichment_results = if (
    length(all_enrichment_results) > 0
  ) {
    combined_enrichments
  } else {
    NULL
  },
  
  module_summary = module_enrichment_summary,
  
  module_status = module_go_status,
  
  errors = go_errors
  
)

saveRDS(
  analysis_summary,
  "./go/go_analysis_summary.rds"
)

# Session information
writeLines(
  capture.output(sessionInfo()),
  "./go/session_info.txt"
)
