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

# Working Directory and Output Setup --------------------------------------

setwd("~/Desktop/Stick_insect/enrichment_analyses")

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
module_assignments <- read_csv(file.path("..", "clitarchus_Module_Summary.csv"))

# Load GO annotations
transID2GO <- read.delim(
  file = "./go/clitarchus_seq_GO.txt",
  sep = "\t",
  header = TRUE,
  fill = TRUE,
  na.strings = "",
  stringsAsFactors = FALSE
) %>%
  dplyr::rename(
    qry_transcript_id = Sequence.Name,
    go_id = Annotation.GO.ID
  ) %>%
  filter(!is.na(go_id), go_id != "")

# Create gene-to-GO mapping list
transID2GOList <- strsplit(as.character(transID2GO$go_id), split = ",")
names(transID2GOList) <- transID2GO$qry_transcript_id

# Gene universe
geneUniverse <- as.character(transID2GO$qry_transcript_id)

# Filter module assignments to genes with GO annotations
module_assignments_go <- module_assignments %>%
  filter(Reference_Transcript_Name %in% geneUniverse)

# Get modules to analyze
modules_to_analyze <- module_assignments_go %>%
  filter(Module_ID != "Module_00") %>%
  group_by(Module_ID) %>%
  summarise(
    genes = list(Reference_Transcript_Name),
    n_genes = n(),
    .groups = 'drop'
  ) %>%
  filter(n_genes >= 10)

# Include Module_00 if sufficient genes
module_00_data <- module_assignments_go %>%
  filter(Module_ID == "Module_00")

if (nrow(module_00_data) >= 10) {
  module_00_summary <- tibble(
    Module_ID = "Module_00",
    genes = list(module_00_data$Reference_Transcript_Name),
    n_genes = nrow(module_00_data)
  )
  modules_to_analyze <- bind_rows(modules_to_analyze, module_00_summary)
}

ontologies <- c("BP", "MF", "CC")

# GO Enrichment Function --------------------------------------------------

perform_module_go_enrichment <- function(module_id, module_genes, ontology = "BP") {
  
  # Create binary gene list
  geneList <- factor(as.integer(geneUniverse %in% module_genes))
  names(geneList) <- geneUniverse
  
  # Check module size
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
    gene2GO = transID2GOList,
    nodeSize = 5
  )
  
  # Run Fisher's exact test
  resultFisher <- runTest(GOdata, algorithm = "parentchild", statistic = "fisher")
  
  # Get GO terms
  allGO <- usedGO(GOdata)
  
  if (length(allGO) == 0) {
    return(NULL)
  }
  
  # Generate results table
  allRes <- GenTable(
    GOdata,
    p_value = resultFisher,
    topNodes = length(allGO),
    numChar = 1000
  )
  
  # Add gene IDs
  allRes$genes <- sapply(allRes$GO.ID, function(go_term) {
    genes_in_term <- genesInTerm(GOdata, go_term)
    module_genes_in_term <- genes_in_term[[1]][genes_in_term[[1]] %in% module_genes]
    paste(module_genes_in_term, collapse = ", ")
  })
  
  # Clean and filter results
  allRes <- allRes %>%
    mutate(
      p_value_clean = case_when(
        p_value == "< 1e-30" ~ 1e-30,
        grepl("^<", p_value) ~ as.numeric(gsub("< ", "", p_value)),
        TRUE ~ as.numeric(p_value)
      ),
      Significant = as.numeric(Significant),
      Expected = as.numeric(Expected),
      module_id = module_id,
      ontology = ontology,
      enrichment_ratio = ifelse(Expected > 0, Significant / Expected, 0)
    ) %>%
    filter(!is.na(p_value_clean), p_value_clean <= 0.05) %>%
    arrange(p_value_clean)
  
  # Replace p_value column
  allRes$p_value <- allRes$p_value_clean
  allRes$p_value_clean <- NULL
  
  return(allRes)
}

# Run GO Enrichment -------------------------------------------------------

all_enrichment_results <- list()

for (ont in ontologies) {
  
  # Collect results for this ontology
  ont_results_list <- list()
  
  for (i in 1:nrow(modules_to_analyze)) {
    module_id <- modules_to_analyze$Module_ID[i]
    module_genes <- modules_to_analyze$genes[[i]]
    
    # Run enrichment analysis
    module_result <- tryCatch({
      perform_module_go_enrichment(module_id, module_genes, ont)
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(module_result) && nrow(module_result) > 0) {
      ont_results_list[[module_id]] <- module_result
    }
  }
  
  # Combine results for this ontology
  if (length(ont_results_list) > 0) {
    ont_results <- bind_rows(ont_results_list)
    all_enrichment_results[[ont]] <- ont_results
    
    # Save individual ontology results
    output_file <- file.path("./go/summary_tables", paste0("All_Modules_", ont, "_Enrichment.csv"))
    write_csv(ont_results, output_file)
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

# Generate Summary Statistics ----------------------------------------------

# Module enrichment summary
module_enrichment_summary <- tibble()

if (length(all_enrichment_results) > 0) {
  module_enrichment_summary <- map_dfr(names(all_enrichment_results), function(ont) {
    ont_data <- all_enrichment_results[[ont]]
    
    if (nrow(ont_data) > 0) {
      ont_data %>%
        group_by(module_id) %>%
        summarise(
          ontology_type = ont,
          n_significant_terms = n(),
          top_p_value = min(p_value, na.rm = TRUE),
          avg_enrichment_ratio = mean(enrichment_ratio, na.rm = TRUE),
          .groups = 'drop'
        )
    } else {
      tibble()
    }
  })
  
  write_csv(
    module_enrichment_summary,
    file.path("./go/summary_tables", "Module_Enrichment_Summary.csv")
  )
} else {
  # Create empty summary file
  empty_summary <- tibble(
    module_id = character(),
    ontology_type = character(),
    n_significant_terms = numeric(),
    top_p_value = numeric(),
    avg_enrichment_ratio = numeric()
  )
  
  write_csv(
    empty_summary,
    file.path("./go/summary_tables", "Module_Enrichment_Summary.csv")
  )
}

# Top enriched terms
if (length(all_enrichment_results) > 0) {
  top_enriched_terms <- bind_rows(all_enrichment_results, .id = "ontology") %>%
    arrange(p_value) %>%
    slice_head(n = 50) %>%
    dplyr::select(ontology, module_id, GO.ID, Term, p_value, enrichment_ratio, Significant)
  
  write_csv(
    top_enriched_terms,
    file.path("./go/summary_tables", "Top_Enriched_Terms_All_Modules.csv")
  )
}

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

# Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./go/session_info.txt")
