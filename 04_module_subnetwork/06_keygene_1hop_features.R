###############################################################
# 005. Key gene 1-hop sub-network feature visualization
###############################################################

rm(list = ls())

library(dplyr)
library(readr)
library(igraph)
library(ggraph)
library(tidyr)
library(stringr)
library(tibble)
library(RColorBrewer)

set.seed(777)

#--------------------------------------------------------------
# 1. Directory
#--------------------------------------------------------------

# === PATHS ===
dir.base <- "/path/to/results"
dir.out  <- file.path(dir.base, "004_Network/module_subnetwork/KeyGene_1hop_networks")
# =============================================
if (!dir.exists(dir.out)) dir.create(dir.out, recursive = TRUE)

#--------------------------------------------------------------
# 2. Module definitions
#--------------------------------------------------------------

source("01_module_definition.R")   # axis_tbl, axis_palette, assign_module()

#--------------------------------------------------------------
# 3. Load data
#--------------------------------------------------------------

df_node <- read_csv(
  file.path(dir.base,
            "004_Network/Leiden_Subgrouping_iterative/<DATE>/Iterative_node_membership_thr200.csv")
)

df_cor <- readRDS(
  file.path(dir.base, "003_Correlation/Filtered_correlation_pairs.rds")
)

# DEG direction matrix (e.g., rows=genes, cols=contrasts, values in {-1,0,1})
file_DEG <- file.path(dir.base, "002_DEG/DEG_matrix_direction.txt")
data_deg <- read.table(file_DEG, check.names = FALSE)

# KeyGene table produced by 002_module_subnetwork.R
node_df <- read_csv(
  file.path(dir.out, "..", "node_centrality_KeyGene.csv")
)

#--------------------------------------------------------------
# 4. Assign module
#--------------------------------------------------------------

df_axis <- assign_module(df_node) %>% filter(!is.na(axis))

#--------------------------------------------------------------
# 4. Build full graph
#--------------------------------------------------------------

th.cor <- 0.9

df_cor_f <- df_cor %>%
  filter(abs(Correlation) > th.cor)

all_axis_genes <- unique(df_axis$node)

edges_axis <- df_cor_f %>%
  filter(Gene1 %in% all_axis_genes,
         Gene2 %in% all_axis_genes)

g_full <- graph_from_data_frame(edges_axis, directed = FALSE)

axis_map <- df_axis %>%
  distinct(node, axis)

V(g_full)$axis <- axis_map$axis[match(V(g_full)$name, axis_map$node)]

#--------------------------------------------------------------
# 5. Add node attributes
#--------------------------------------------------------------

# KeyGene information
V(g_full)$KeyGene <- node_df$KeyGene[match(V(g_full)$name, node_df$name)]
V(g_full)$KeyGene[is.na(V(g_full)$KeyGene)] <- FALSE

# Prepare DEG data as dataframe for easier joining
deg_df <- data_deg %>%
  rownames_to_column("gene") %>%
  as_tibble()

deg_cols <- colnames(data_deg)

# DEG column name mapping
deg_col_mapping <- c(
  "T25_H1_T42_H1" = "25℃ (1h) vs. 42℃ (1h)",
  "T25_D1_T42_D1" = "25℃ (1d) vs. 42℃ (1d)",
  "T25_H1_T25_D1" = "25℃ (1h) vs. 25℃ (1d)",
  "T42_H1_T42_D1" = "42℃ (1h) vs. 42℃ (1d)"
)

deg_col_mapping <- c(
  "T25_H1_T42_H1" = "(1h) 25℃ vs. 42℃",
  "T25_D1_T42_D1" = "(1d) 25℃ vs. 42℃",
  "T25_H1_T25_D1" = "(25℃) 1h vs. 1d",
  "T42_H1_T42_D1" = "(42℃) 1h vs. 1d"
)


#--------------------------------------------------------------
# 6. Get Key Genes
#--------------------------------------------------------------

key_genes <- node_df %>%
  filter(KeyGene == TRUE) %>%
  pull(name)

cat("Total Key Genes:", length(key_genes), "\n")

#--------------------------------------------------------------
# 7. Color palettes
#--------------------------------------------------------------

pal_axis <- axis_palette   # from 000_module_definition.R

pal_deg <- c(
  "-1" = "#4575B4",  # Down (blue)
  "1"  = "#D73027"   # Up (red)
)

#--------------------------------------------------------------
# 8. Function: Extract 1-hop subgraph and save info
#--------------------------------------------------------------

extract_1hop_and_save <- function(graph, gene_name, axis_map, deg_df, output_dir) {
  if (!gene_name %in% V(graph)$name) {
    return(NULL)
  }
  
  # Get neighbors
  neighbors <- neighbors(graph, gene_name)
  
  # Include focal gene + neighbors
  nodes_to_keep <- c(gene_name, neighbors$name)
  
  # Extract subgraph
  subg <- induced_subgraph(graph, nodes_to_keep)
  
  # Mark focal gene
  V(subg)$is_focal <- V(subg)$name == gene_name
  
  # Create info table
  info_df <- tibble(
    gene = V(subg)$name,
    is_focal = V(subg)$is_focal,
    axis = V(subg)$axis,
    KeyGene = V(subg)$KeyGene
  ) %>%
    left_join(deg_df, by = "gene")
  
  # Save info table
  write_csv(info_df, file.path(output_dir, paste0(gene_name, "_1hop_info.csv")))
  
  return(list(subg = subg, info = info_df))
}

#--------------------------------------------------------------
# 9. Function: Plot 1-hop network (Cluster info)
#--------------------------------------------------------------

plot_1hop_cluster <- function(subg, gene_name) {
  if (is.null(subg) || vcount(subg) == 0) {
    return(NULL)
  }
  
  # Count nodes
  n_nodes <- vcount(subg)
  
  # Convert KeyGene to factor with desired order
  V(subg)$KeyGene_factor <- factor(V(subg)$KeyGene, levels = c(TRUE, FALSE))
  
  set.seed(123)
  
  p <- ggraph(subg, layout = "fr") +
    geom_edge_link(alpha = 0.3, color = "grey50", width = 0.5) +
    geom_node_point(
      aes(fill = axis, 
          shape = KeyGene_factor,
          size = is_focal),
      color = "black", 
      stroke = 0.5, 
      alpha = 0.8
    ) +
    geom_node_text(
      aes(label = ifelse(is_focal, name, "")),
      size = 4, 
      fontface = "bold",
      vjust = -1.5
    ) +
    scale_fill_manual(
      values = pal_axis,
      name = "Cluster",
      guide = guide_legend(override.aes = list(size = 5, shape = 21))
    ) +
    scale_shape_manual(
      values = c("TRUE" = 24, "FALSE" = 21),
      name   = "Gene type",
      labels = c("KeyGene", "Other"),
      guide = guide_legend(override.aes = list(size = 5))
    ) +
    scale_size_manual(
      values = c("FALSE" = 3, "TRUE" = 6),
      guide = "none"
    ) +
    guides(
	  shape = guide_legend(order = 1, override.aes = list(size = 5)),
      fill = guide_legend(order = 2, override.aes = list(size = 5, shape = 21, color = "black", stroke = 0.5))
    ) +
    theme_graph() +
    labs(title = paste0(gene_name, " - 1-hop Network (Cluster) [n=", n_nodes, "]")) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9)
    )
  
  return(p)
}

#--------------------------------------------------------------
# 10. Function: Plot 1-hop network (DEG info)
#--------------------------------------------------------------

plot_1hop_deg <- function(info_df, gene_name, deg_col, subg, deg_col_mapping) {
  if (is.null(subg) || vcount(subg) == 0) {
    return(NULL)
  }
  
  # Count nodes
  n_nodes <- vcount(subg)
  
  # Get DEG values from info_df
  deg_values <- info_df[[deg_col]]
  
  # Convert 0 to NA (no change -> NA)
  deg_values[deg_values == 0] <- NA
  
  # Add to subgraph
  V(subg)$deg_value <- deg_values
  # Change order: Up (1), Down (-1), NA
  V(subg)$deg_factor <- factor(deg_values, levels = c(1, -1))
  
  # Convert KeyGene to factor with desired order
  V(subg)$KeyGene_factor <- factor(V(subg)$KeyGene, levels = c(TRUE, FALSE))
  
  # Get subtitle
  subtitle_text <- deg_col_mapping[deg_col]
  
  set.seed(123)
  
  p <- ggraph(subg, layout = "fr") +
    geom_edge_link(alpha = 0.3, color = "grey50", width = 0.5) +
    geom_node_point(
      aes(fill = deg_factor,
          shape = KeyGene_factor,
          size = is_focal),
      color = "black", 
      stroke = 0.5, 
      alpha = 0.8
    ) +
    geom_node_text(
      aes(label = ifelse(is_focal, name, "")),
      size = 4, 
      fontface = "bold",
      vjust = -1.5
    ) +
    scale_fill_manual(
      values = c("1" = "#D73027", "-1" = "#4575B4"),
      name = "DEG",
      labels = c("1" = "Up", "-1" = "Down"),
      na.value = "white"
    ) +
    scale_shape_manual(
      values = c("TRUE" = 24, "FALSE" = 21),
      name   = "Gene type",
      labels = c("KeyGene", "Other")
    ) +
    scale_size_manual(
      values = c("FALSE" = 3, "TRUE" = 6),
      guide = "none"
    ) +
    guides(
      shape = guide_legend(order = 1, override.aes = list(size = 5)),
      fill = guide_legend(order = 2, override.aes = list(size = 5, shape = 21, color = "black", stroke = 0.5))
    ) +
    theme_graph() +
    labs(
      title = paste0(gene_name, " - 1-hop Network [n=", n_nodes, "]"),
      subtitle = subtitle_text
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9)
    )
  
  return(p)
}

#--------------------------------------------------------------
# 11. Generate plots for all key genes
#--------------------------------------------------------------

deg_comparisons <- colnames(data_deg)

cat("\nGenerating 1-hop networks for", length(key_genes), "key genes...\n")
cat("Each gene will have 5 plots + 1 info CSV\n\n")

# Progress tracking
pb <- txtProgressBar(min = 0, max = length(key_genes), style = 3)

success_count <- 0
error_count <- 0


for (i in seq_along(key_genes)) {
  gene <- key_genes[i]
  
  tryCatch({
    # Create output directory for this gene
    gene_dir <- file.path(dir.out, gene)
    if (!dir.exists(gene_dir)) dir.create(gene_dir, recursive = TRUE)
    
    # Extract 1-hop subgraph and save info
    result <- extract_1hop_and_save(g_full, gene, axis_map, deg_df, gene_dir)
    
    if (is.null(result)) {
      cat("\nWarning: Could not extract subgraph for", gene, "\n")
      error_count <- error_count + 1
      next
    }
    
    subg <- result$subg
    info_df <- result$info
    
    # Plot 1: Cluster information
    p_cluster <- plot_1hop_cluster(subg, gene)
    if (!is.null(p_cluster)) {
      ggsave(
        file.path(gene_dir, paste0(gene, "_cluster.png")),
        p_cluster, 
        width = 8, 
        height = 7, 
        dpi = 300
      )
    }
    
    # Plots 2-5: DEG information
    for (deg_col in deg_comparisons) {
      p_deg <- plot_1hop_deg(info_df, gene, deg_col, subg, deg_col_mapping)
      if (!is.null(p_deg)) {
        deg_col_clean <- gsub("_", "-", deg_col)
        ggsave(
          file.path(gene_dir, paste0(gene, "_", deg_col_clean, ".png")),
          p_deg, 
          width = 8, 
          height = 7, 
          dpi = 300
        )
      }
    }
    
    success_count <- success_count + 1
    
  }, error = function(e) {
    cat("\nError processing", gene, ":", e$message, "\n")
    error_count <- error_count + 1
  })
  
  # Update progress
  setTxtProgressBar(pb, i)
}

close(pb)

#--------------------------------------------------------------
# 12. Summary statistics
#--------------------------------------------------------------

cat("\n\n===== Summary =====\n")
cat("Total key genes:", length(key_genes), "\n")
cat("Successfully processed:", success_count, "\n")
cat("Errors:", error_count, "\n")
cat("Output directory:", dir.out, "\n")
cat("Files per gene: 6 (5 plots + 1 info CSV)\n")
cat("Total files generated:", success_count * 6, "\n")
cat("\n===== ALL DONE =====\n")