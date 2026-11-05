################################################################################
#
# Title: WGCNA Module KEGG Enrichment Analysis
#        Colorado Potato Beetle (CPB) Diapause Gene Co-expression Networks
#
# Author: Alex Torson
# Institution: USDA-ARS
# Email: Alex.Torson@usda.gov
#
# Description: KEGG pathway enrichment analysis on WGCNA modules to identify
#              biological pathways significantly enriched in each co-expression
#              module, providing functional interpretation of network structure.
#
# Dependencies: Must run 01_core_wgcna_analysis.R first
#
# Input files:
#   - 01_core_results/Module_Summary.csv
#   - Annotation/cpb_diapause_KAAS.txt
#   - Annotation/cpb_diapause_annotation.csv
#
# Output files:
#   - 13_kegg_enrichment/module_enrichments/*.csv (individual module results)
#   - 13_kegg_enrichment/summary_tables/*.csv (combined results)
#   - 13_kegg_enrichment/plots/*.png (visualization plots)
#
################################################################################

# Package Loading and Setup -----------------------------------------------

library(tidyverse)
library(pathview)
library(gage)
library(clusterProfiler)

httr::set_config(httr::config(ssl_verifypeer = FALSE))

# Working Directory and Output Setup --------------------------------------

setwd("~/Desktop/Stick_insect/Stick_insect_retry/enrichment_analyses")

kegg_dirs <- c(
  "kegg/module_enrichments",
  "kegg/summary_tables", 
  "kegg/plots",
  "kegg/pathway_visualizations"
)

walk(kegg_dirs, ~ if (!dir.exists(.x)) dir.create(.x, recursive = TRUE))

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
    legend.position = "bottom"
  )

# Data Loading and Preparation --------------------------------------------

# Load WGCNA module assignments from script 01
module_assignments <- read_csv(file.path("~/Desktop/Stick_insect/clitarchus_Module_Summary_Annotated.csv"))

# Load BLAST annotations from BRR species
blast_annot <- read.csv(file.path("~/Desktop/Stick_insect/clitarchus_annotation.csv"))

# Load KEGG ortholog mappings from KAAS annotation
cli_kegg <- read.table(
  "~/Desktop/Stick_insect/clitarchus_KAAS_query.txt",
  header = FALSE,
  colClasses = c("character", "factor"),
  sep = "",
  fill = TRUE,
  col.names = c("qry_transcript_id", "kegg_id"),
  na.strings = ""
) %>%
  drop_na()

# Load gene annotations
library(rtracklayer)
cli_annot <- import("~/Desktop/Stick_insect/cli_annevo_gemoma_combined.gtf")
cli_annot <- as.data.frame(cli_annot) # convert to data frame

# Merge KEGG IDs with gene annotations
cli_annot <- cli_annot %>%
  left_join(
    cli_kegg,
    by = c("transcript_id" = "qry_transcript_id")
  ) %>%
  dplyr::rename(qry_transcript_id = transcript_id)

# Merge ref_transcript_name from BLAST to cli_annot
cli_annot <- cli_annot %>%
  left_join(
    blast_annot %>% select(qry_gene_id, gene_annotation),
    by = c("gene_id" = "qry_gene_id")
  )

# Create gene-to-KEGG ortholog mapping
kegg_genes <- cli_annot %>%
  dplyr::select(gene_id, gene_annotation, kegg_id) %>%
  drop_na(gene_id, kegg_id) %>%
  distinct()

# -------------------------------------------------------
# DEFINE EXPRESSED GENE UNIVERSE
# -------------------------------------------------------

# Include all genes from the original WGCNA network:
# 18 retained modules + excluded modules + grey/unassigned

expressed_gene_ids <- module_assignments %>%
  pull(gene_id) %>%
  unique()

# Confirm expected number of expressed genes
length(expressed_gene_ids)

stopifnot(length(expressed_gene_ids) == 19252)

# Filter KEGG annotations to expressed genes only
kegg_genes <- kegg_genes %>%
  filter(gene_id %in% expressed_gene_ids) %>%
  distinct(gene_id, kegg_id, .keep_all = TRUE)

# Define background universe using unique KO IDs
universe_kegg <- kegg_genes %>%
  pull(kegg_id) %>%
  unique() %>%
  as.character()

# Download KEGG pathway-to-KO mappings
kegg_data <- clusterProfiler::download_KEGG(
  species = "ko",
  keggType = "KEGG",
  keyType = "kegg"
)

# Create pathway-to-KO mapping
kegg.gs <- split(
  kegg_data$KEGGPATHID2EXTID$to,
  kegg_data$KEGGPATHID2EXTID$from
)

# Create pathway ID-to-name mapping
kegg.gs_names <- kegg_data$KEGGPATHID2NAME

# Check objects
length(kegg.gs)
head(kegg.gs)
head(kegg.gs_names)

# Filter module assignments to genes with KEGG annotations
module_assignments_kegg <- module_assignments %>%
  filter(gene_id %in% kegg_genes$gene_id)

# Check numbers
cat("Total expressed genes:", length(expressed_gene_ids), "\n")
cat("Expressed genes with KEGG annotations:",
    n_distinct(kegg_genes$gene_id), "\n")
cat("Unique KO IDs in background universe:",
    length(universe_kegg), "\n")

# KEGG Enrichment Functions -----------------------------------------------

# Test KEGG pathway enrichment for a single module
perform_module_kegg_enrichment <- function(module_id, module_genes) {
  
  # Get KEGG IDs for genes in this module
  module_kegg <- tibble(gene_id = module_genes) %>%
    left_join(kegg_genes %>% select(gene_id, kegg_id), by = "gene_id") %>%
    drop_na(kegg_id) %>%
    distinct()
  
  # Extract unique KEGG ortholog IDs - contained in expressed gene set
  module_kegg_ids <- intersect(
    unique(as.character(module_kegg$kegg_id)),
    universe_kegg
  )
  
  # Skip modules with fewer than 5 annotated KO IDs
  if (length(module_kegg_ids) < 5) {
    return(NULL)
  }
  
  # Over-representation analysis using Fisher's exact test
  enrich_result <- enrichKEGG(
    gene = module_kegg_ids,
    universe = universe_kegg,
    organism = "ko",
    keyType = 'kegg',
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
  
  # Check whether enrichment returned any results
  if (is.null(enrich_result) || nrow(enrich_result@result) == 0) {
    return(NULL)
  }
  
  # Convert to dataframe and and retain FDR-significant pathways
  enrich_df <- as.data.frame(enrich_result) %>%
    filter(p.adjust <= 0.05) %>%
    mutate(
      module_id = module_id,
      bg_ratio_numeric =
        as.numeric(sub("/.*", "", BgRatio)) /
        as.numeric(sub(".*/", "", BgRatio)),
      enrichment_ratio =
        (Count / length(module_kegg_ids)) / bg_ratio_numeric
    ) %>%
    arrange(p.adjust)
  
  return(enrich_df)
}

# Generate pathway-gene mappings for a module
generate_module_pathway_mappings <- function(module_id, module_genes) {
  
  # Get KEGG annotations for module genes
  module_kegg <- tibble(gene_id = module_genes) %>%
    left_join(
      kegg_genes,
      by = "gene_id"
    ) %>%
    drop_na(kegg_id) %>%
    distinct()
  
  if (nrow(module_kegg) == 0) {
    return(NULL)
  }
  
  # Create pathway-to-KO mapping
  pathway_kegg <- tibble(
    pathway_id = rep(names(kegg.gs), lengths(kegg.gs)),
    kegg_id = unlist(kegg.gs, use.names = FALSE)
  )
  
  # Add pathway names
  pathway_kegg <- pathway_kegg %>%
    mutate(
      pathway_name = kegg.gs_names$to[
        match(pathway_id, kegg.gs_names$from)
      ]
    )
  
  # Match module genes to KEGG pathways
  pathway_transcripts <- pathway_kegg %>%
    inner_join(
      module_kegg,
      by = "kegg_id",
      relationship = "many-to-many"
    ) %>%
    select(
      pathway_id,
      pathway_name,
      gene_id,
      gene_annotation,
      kegg_id
    ) %>%
    distinct()
  
  if (nrow(pathway_transcripts) == 0) {
    return(NULL)
  }
  
  # Create pathway gene count summary
  pathway_gene_counts <- pathway_transcripts %>%
    group_by(pathway_id, pathway_name) %>%
    summarise(
      gene_count = n_distinct(gene_id),
      KO_count = n_distinct(kegg_id),
      .groups = "drop"
    ) %>%
    arrange(desc(gene_count))
  
  return(list(
    counts = pathway_gene_counts,
    mappings = pathway_transcripts
  ))
}

# Run KEGG Enrichment Analysis --------------------------------------------

# Get modules to analyze (require minimum size)
modules_to_analyze <- module_assignments_kegg %>%
  filter(Module_Type == "retained") %>%
  group_by(Module_ID) %>%
  summarise(
    genes = list(gene_id),
    n_genes = n(),
    .groups = 'drop'
  ) %>%
  filter(n_genes >= 10)

# Run enrichment analysis for all modules
all_kegg_enrichments <- list()
all_pathway_mappings <- list()

for (i in 1:nrow(modules_to_analyze)) {
  module_id <- modules_to_analyze$Module_ID[i]
  module_genes <- modules_to_analyze$genes[[i]]
  
  # KEGG pathway enrichment
  enrichment_result <- perform_module_kegg_enrichment(module_id, module_genes)
  if (!is.null(enrichment_result) && nrow(enrichment_result) > 0) {
    all_kegg_enrichments[[module_id]] <- enrichment_result
  }
  
  # Pathway gene mappings
  pathway_result <- generate_module_pathway_mappings(module_id, module_genes)
  if (!is.null(pathway_result)) {
    all_pathway_mappings[[module_id]] <- pathway_result
  }
}


# Save Individual Module Results ------------------------------------------

for (i in seq_len(nrow(modules_to_analyze))) {
  
  module_id <- modules_to_analyze$Module_ID[i]
  module_genes <- modules_to_analyze$genes[[i]]
  
  # Create module directory
  module_dir <- file.path("./kegg/module_enrichments", module_id)
  dir.create(module_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save enrichment results
  enrich_df <- all_kegg_enrichments[[module_id]]
  
  if (!is.null(enrich_df) && nrow(enrich_df) > 0) {
    
    write_csv(
      enrich_df,
      file.path(module_dir, paste0(module_id, "_KEGG_enrichment.csv"))
    )
    
    # Create pathway-gene mapping from enriched pathways
    pathway_gene_mapping <- map_dfr(seq_len(nrow(enrich_df)), function(j) {
      
      pathway_id <- enrich_df$ID[j]
      pathway_name <- enrich_df$Description[j]
      
      # Extract KO IDs from enriched pathway
      pathway_kos <- strsplit(enrich_df$geneID[j], "/", fixed = TRUE)[[1]]
      
      # Map KO IDs back to original genes
      tibble(kegg_id = pathway_kos) %>%
        inner_join(
          kegg_genes %>%
            select(gene_id, gene_annotation, kegg_id) %>%
            distinct(),
          by = "kegg_id"
        ) %>%
        mutate(
          pathway_id = pathway_id,
          pathway_name = pathway_name
        ) %>%
        select(pathway_id, pathway_name, gene_id, gene_annotation, kegg_id)
      
    }) %>%
      distinct() %>%
      arrange(pathway_id, gene_id)
    
    # Save pathway-gene mapping
    write_csv(
      pathway_gene_mapping,
      file.path(module_dir, paste0(module_id, "_pathway_gene_ids.csv"))
    )
    
  } else {
    
    write_lines(
      paste("No significant KEGG pathways found for", module_id),
      file.path(module_dir, paste0(module_id, "_no_KEGG_enrichment.txt"))
    )
  }
  
  # Save pathway mappings if available
  pathway_mapping <- all_pathway_mappings[[module_id]]
  
  if (!is.null(pathway_mapping)) {
    
    write_csv(
      pathway_mapping$counts,
      file.path(module_dir, paste0(module_id, "_pathway_gene_counts.csv"))
    )
    
    write_csv(
      pathway_mapping$mappings,
      file.path(module_dir, paste0(module_id, "_pathway_gene_mappings.csv"))
    )
  }
  
  # Save module gene list
  write_lines(
    module_genes,
    file.path(module_dir, paste0(module_id, "_gene_list.txt"))
  )
}


# Create Summary Tables ---------------------------------------------------

# Combined enrichment results across all modules
if (length(all_kegg_enrichments) > 0) {
  
  combined_enrichments <- bind_rows(
    all_kegg_enrichments,
    .id = "module_id"
  )
  
  if (nrow(combined_enrichments) > 0) {
    
    write_csv(
      combined_enrichments,
      file.path("./kegg/summary_tables", "All_Modules_KEGG_Enrichment.csv")
    )
    
    # Module enrichment summary
    module_kegg_summary <- combined_enrichments %>%
      group_by(module_id) %>%
      summarise(
        n_significant_pathways = n(),
        top_p_adjust = min(p.adjust, na.rm = TRUE),
        mean_enrichment_ratio = mean(enrichment_ratio, na.rm = TRUE),
        top_pathway = Description[which.min(p.adjust)],
        .groups = "drop"
      ) %>%
      arrange(desc(n_significant_pathways))
    
    write_csv(
      module_kegg_summary,
      file.path("./kegg/summary_tables", "Module_KEGG_Summary.csv")
    )
    
    # Top enriched pathways across all modules
    top_pathways <- combined_enrichments %>%
      arrange(p.adjust) %>%
      slice_head(n = 30) %>%
      select(
        module_id, Description, pvalue, p.adjust,
        qvalue, Count, enrichment_ratio
      )
    
    write_csv(
      top_pathways,
      file.path("./kegg/summary_tables", "Top_KEGG_Pathways_All_Modules.csv")
    )
  }
}


# Combined Pathway Mappings ------------------------------------------------

if (length(all_pathway_mappings) > 0) {
  
  # Combine pathway gene counts
  combined_pathway_counts <- map_dfr(
    all_pathway_mappings,
    ~ .x$counts,
    .id = "module_id"
  )
  
  # Combine detailed pathway-gene mappings
  combined_pathway_mappings <- map_dfr(
    all_pathway_mappings,
    ~ .x$mappings,
    .id = "module_id"
  )
  
  write_csv(
    combined_pathway_counts,
    file.path("./kegg/summary_tables", "All_Modules_Pathway_Counts.csv")
  )
  
  write_csv(
    combined_pathway_mappings,
    file.path("./kegg/summary_tables", "All_Modules_Pathway_Mappings.csv")
  )
}


# Create Visualization Plots ----------------------------------------------

if (length(all_kegg_enrichments) > 0) {
  
  # Module enrichment heatmap
  if (exists("module_kegg_summary") && nrow(module_kegg_summary) > 0) {
    
    p_module_heatmap <- ggplot(module_kegg_summary, 
                               aes(x = 1, y = module_id, fill = n_significant_pathways)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = n_significant_pathways), color = "black", size = 3) +
      scale_fill_gradient(low = "white", high = sig_colors[["TRUE"]], name = "Significant\nPathways") +
      labs(
        x = "",
        y = "WGCNA Module"
      ) +
      theme_cpb +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right"
      )
    
    ggsave(
      file.path("./kegg/plots", "Module_KEGG_Enrichment_Heatmap.pdf"),
      p_module_heatmap, width = 6, height = 8, dpi = 300
    )
    
    ggsave(
      file.path("./kegg/plots", "Module_KEGG_Enrichment_Heatmap.png"),
      p_module_heatmap, width = 6, height = 8, dpi = 300
    )
  }
  
  # Top pathways bar plot
  if (exists("top_pathways") && nrow(top_pathways) > 0) {
    
    top_pathways_plot <- top_pathways %>%
      slice_head(n = 20) %>%
      mutate(
        pathway_label = paste0(str_trunc(Description, 40), " (", module_id, ")"),
        neg_log_p = -log10(pvalue)
      )
    
    p_top_pathways <- ggplot(top_pathways_plot, 
                             aes(x = reorder(pathway_label, neg_log_p), y = neg_log_p)) +
      geom_col(fill = sig_colors[["TRUE"]], alpha = 0.7) +
      coord_flip() +
      labs(
        x = "KEGG Pathway (Module)",
        y = "-log10(p-value)"
      ) +
      theme_cpb +
      theme(
        axis.text.y = element_text(size = 8),
        legend.position = "none"
      )
    
    ggsave(
      file.path("./kegg/plots", "Top_KEGG_Pathways.pdf"),
      p_top_pathways, width = 12, height = 8, dpi = 300
    )
    
    ggsave(
      file.path("./kegg/plots", "Top_KEGG_Pathways.png"),
      p_top_pathways, width = 12, height = 8, dpi = 300
    )
  }
  
  # Enrichment ratio vs significance scatter plot
  if (exists("combined_enrichments") && nrow(combined_enrichments) > 0) {
    
    plot_data <- combined_enrichments %>%
      filter(
        !is.na(enrichment_ratio), 
        !is.infinite(enrichment_ratio),
        enrichment_ratio > 0,
        !is.na(pvalue),
        pvalue > 0
      )
    
    if (nrow(plot_data) > 0) {
      p_enrichment_scatter <- ggplot(plot_data, 
                                     aes(x = enrichment_ratio, y = -log10(pvalue), color = module_id)) +
        geom_point(alpha = 0.7, size = 2) +
        labs(
          x = "Enrichment Ratio",
          y = "-log10(p-value)",
          color = "Module"
        ) +
        theme_cpb +
        theme(legend.position = "right")
      
      ggsave(
        file.path("./kegg/plots", "Enrichment_Ratio_vs_Significance.pdf"),
        p_enrichment_scatter, width = 10, height = 6, dpi = 300
      )
      
      ggsave(
        file.path("./kegg/plots", "Enrichment_Ratio_vs_Significance.png"),
        p_enrichment_scatter, width = 10, height = 6, dpi = 300
      )
    }
  }
}

# Generate Pathway Visualizations (Top Pathways Only) --------------------

if (length(all_kegg_enrichments) > 0) {
  
  # Get top 5 most significant pathways for visualization
  # Rara: changed arrange(pvalue) to arrange(p.adjust)
  top_pathways_for_viz <- bind_rows(all_kegg_enrichments, .id = "module_id") %>%
    arrange(p.adjust) %>%
    slice_head(n = 5) %>%
    dplyr::select(module_id, ID, Description)
  
  # Exclude problematic pathways that don't visualize well
#  excluded_pathways <- c("ko04723")
#    "ko01100", "ko01110", "ko04723", "ko01120", "ko04215",
#    "ko00510", "ko00511", "ko00512", "ko00513", "ko00514", "ko00515",
#    "ko00531", "ko00532", "ko00533", "ko00534", "ko00563",
#    "ko00601", "ko00603", "ko00604", "ko00966",
#    "ko01200", "ko01040", "ko01210", "ko01212", "ko01220",
#    "ko01230", "ko01232", "ko01240", "ko01250"
#  )
  
#  top_pathways_for_viz <- top_pathways_for_viz %>%
#    filter(!ID %in% excluded_pathways)
  
  if (nrow(top_pathways_for_viz) > 0) {
    
    viz_dir <- file.path("./kegg/pathway_visualizations")
    
    for (i in 1:nrow(top_pathways_for_viz)) {
      module_id <- top_pathways_for_viz$module_id[i]
      pathway_id <- top_pathways_for_viz$ID[i]
      pathway_desc <- top_pathways_for_viz$Description[i]
      
      # Get module genes with KEGG annotations
      module_genes <- modules_to_analyze$genes[[which(modules_to_analyze$Module_ID == module_id)]]
      
      module_kegg <- tibble(gene_id = module_genes) %>%
        left_join(kegg_genes, by = "gene_id") %>%
        drop_na(kegg_id) %>%
        distinct(kegg_id, .keep_all = TRUE)
      
      if (nrow(module_kegg) == 0) next
      
      # Create binary presence vector for module genes
      module_kegg_presence <- setNames(rep(1, nrow(module_kegg)), module_kegg$kegg_id)
      
      # Create pathway visualization
      setwd(viz_dir)
      
      pathway_clean_name <- gsub("[^a-zA-Z0-9_]", "_", pathway_desc)
      out_suffix <- paste0(module_id, "_", pathway_clean_name)
      
      # Rara: removed discrete = list(gene = TRUE) AND limit = list(gene = 1)
      pathview(
        gene.data = module_kegg_presence,
        pathway.id = pathway_id,
        species = "ko",
        out.suffix = out_suffix,
        kegg.native = TRUE,
        low = "white",
        high = "red",
        same.layer = FALSE
      )
      
      # Return to working directory
      setwd("~/Desktop/Stick_insect/Stick_insect_retry/enrichment_analyses")
    }
  }
}

# Save Analysis Summary Objects -------------------------------------------

analysis_summary <- list(
  modules_analyzed = nrow(modules_to_analyze),
  modules_with_enrichments = length(all_kegg_enrichments),
  total_genes_with_kegg = n_distinct(kegg_genes$gene_id),
  genes_in_modules_with_kegg = nrow(module_assignments_kegg),
  enrichment_results = if (length(all_kegg_enrichments) > 0) combined_enrichments else NULL,
  pathway_mappings = if (length(all_pathway_mappings) > 0) combined_pathway_mappings else NULL
)

saveRDS(analysis_summary, file.path("./kegg", "kegg_analysis_summary.rds"))

# Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./kegg/session_info.txt")

