###############################################################
# Full-network non-DEG KeyGene cross-module neighbor analysis
###############################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(openxlsx)
})

## ------------------------------------------------------------------
## 0. User settings
## ------------------------------------------------------------------

# Put genes here when you want to inspect specific genes.
# Example: TARGET_GENES <- c("g35527", "g813")
TARGET_GENES <- character(0)

# TRUE: if TARGET_GENES is empty, plot every full-network KeyGene that is
# not in the DEG union. FALSE: use the top N ranked genes only.
PROCESS_ALL_NONDEG_KEYGENES <- TRUE

# Used only when PROCESS_ALL_NONDEG_KEYGENES is FALSE.
N_AUTO_GENES <- 2

# Four-image set per gene:
#   1) module, 2) heat 1h, 3) heat 1d, 4) 42C time comparison.
PLOT_COMPS <- c("T25_H1_T42_H1", "T25_D1_T42_D1", "T42_H1_T42_D1")

# When the selected gene has too many neighbors, thresholds are tried
# in this order until fewer than MAX_NEIGHBORS connected genes remain.
COR_THRESHOLDS <- c(0.90, 0.95, 0.99)
MAX_NEIGHBORS <- 100

# === EDIT THESE PATHS FOR YOUR ENVIRONMENT ===
paths <- list(
  full_base    = "/path/to/full_network/3_Results",
  degonly_base = "/path/to/Network_DEGonly",
  output_base  = "/path/to/output_root"
)
# =============================================

paths$node_membership <- file.path(
  paths$full_base,
  "004_Network/Leiden_Subgrouping_iterative/2025-10-28/Iterative_node_membership_thr200.csv"
)
paths$cor_pairs <- file.path(
  paths$full_base,
  "003_Correlation/Filtered_correlation_pairs_30samples.rds"
)
paths$full_keygene <- file.path(
  paths$full_base,
  "004_Network/20251213_further_subcluster_heat_photo2/node_centrality_KeyGene_max_membership.csv"
)
paths$deg_rdata <- file.path(
  paths$degonly_base,
  "Nicotiana_DEGs_Heat_2025-11-12_signicant.RData"
)
paths$tpm <- file.path(paths$degonly_base, "Nb_ymkim_TPM.txt")

out_dir <- file.path(
  paths$output_base,
  paste0("Full_nonDEG_KeyGene_balanced_modules_", format(Sys.time(), "%Y-%m-%d_%H%M%S"))
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", out_dir)

## ------------------------------------------------------------------
## 1. Helpers
## ------------------------------------------------------------------

sanitize_df <- function(df) {
  df %>%
    mutate(across(where(is.character), ~ iconv(.x, from = "", to = "UTF-8", sub = "")))
}

safe_div <- function(num, den) {
  ifelse(is.na(den) | den == 0, NA_real_, num / den)
}

read_deg_union_from_rdata <- function(rdata_path) {
  deg_env <- new.env(parent = emptyenv())
  load(rdata_path, envir = deg_env)
  if (!exists("result_list", envir = deg_env)) {
    stop("result_list was not found in: ", rdata_path)
  }

  result_list <- get("result_list", envir = deg_env)
  target_comps <- names(result_list)[1:2]

  deg_long <- bind_rows(lapply(target_comps, function(comp) {
    as.data.frame(result_list[[comp]]) %>%
      rownames_to_column("gene") %>%
      filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 2) %>%
      mutate(
        comparison = comp,
        DEG_direction = case_when(
          log2FoldChange > 0 ~ "Up",
          log2FoldChange < 0 ~ "Down",
          TRUE ~ NA_character_
        )
      )
  }))

  plot_comps <- intersect(PLOT_COMPS, names(result_list))

  direction_wide <- bind_rows(lapply(plot_comps, function(comp) {
    as.data.frame(result_list[[comp]]) %>%
      rownames_to_column("gene") %>%
      transmute(
        gene,
        comparison = comp,
        DEG_value = case_when(
          !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 2 & log2FoldChange > 0 ~ 1,
          !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 2 & log2FoldChange < 0 ~ -1,
          TRUE ~ 0
        ),
        log2FoldChange = log2FoldChange,
        padj = padj
      )
  })) %>%
    select(gene, comparison, DEG_value) %>%
    pivot_wider(names_from = comparison, values_from = DEG_value, values_fill = 0)

  list(
    target_comps = target_comps,
    plot_comps = plot_comps,
    deg_long = deg_long,
    deg_union = unique(deg_long$gene),
    direction_wide = direction_wide
  )
}

make_axis_map <- function(node_membership_path) {
  axis_tbl <- tribble(
    ~axis_code, ~pattern,
    "PS",   "^C_1_4_|^C_1_1_|^C_6_1_1",
    "HEAT", "^C_6_2_|^C_6_3|^C_6_4|^C_1_2_2_",
    "ROS",  "^C_1_3_1_1_1|^C_1_3_2_"
  )

  read_csv(node_membership_path, show_col_types = FALSE) %>%
    mutate(
      axis = case_when(
        str_detect(subcluster_id, axis_tbl$pattern[axis_tbl$axis_code == "PS"]) ~ "PS",
        str_detect(subcluster_id, axis_tbl$pattern[axis_tbl$axis_code == "HEAT"]) ~ "HEAT",
        str_detect(subcluster_id, axis_tbl$pattern[axis_tbl$axis_code == "ROS"]) ~ "ROS",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(axis)) %>%
    distinct(node, axis, subcluster_id)
}

ensure_module_cols <- function(df, module_cols = c("HEAT", "ROS", "PS")) {
  for (module in module_cols) {
    if (!module %in% colnames(df)) df[[module]] <- 0L
  }
  df
}

collapse_genes <- function(x) {
  x <- sort(unique(x[!is.na(x)]))
  if (!length(x)) return("")
  paste(x, collapse = ";")
}

get_gene_edges_at_threshold <- function(cor_dt, gene, threshold) {
  edges1 <- cor_dt[Gene1 == gene & abs_cor >= threshold,
                   .(from = Gene1, to = Gene2, Correlation, abs_cor)]
  edges2 <- cor_dt[Gene2 == gene & abs_cor >= threshold,
                   .(from = Gene2, to = Gene1, Correlation, abs_cor)]
  rbindlist(list(edges1, edges2), use.names = TRUE)
}

choose_threshold <- function(cor_dt, gene, thresholds, max_neighbors) {
  counts <- sapply(thresholds, function(th) {
    uniqueN(get_gene_edges_at_threshold(cor_dt, gene, th)$to)
  })

  hit <- which(counts > 0 & counts < max_neighbors)
  if (length(hit)) {
    chosen_idx <- hit[1]
  } else {
    positive <- which(counts > 0)
    chosen_idx <- if (length(positive)) positive[length(positive)] else length(thresholds)
  }

  tibble(
    threshold = thresholds[chosen_idx],
    n_neighbors = counts[chosen_idx],
    tried_thresholds = paste0(thresholds, "=", counts, collapse = "; ")
  )
}

plot_module_network <- function(layout_df, graph, outfile, gene, threshold) {
  pal_axis <- c(
    "PS" = "#7FBF7F",
    "HEAT" = "#FB9A99",
    "ROS" = "#EDD072"
  )

  p <- ggraph(layout_df) +
    geom_edge_link(aes(width = abs_cor), alpha = 0.28, color = "grey45") +
    geom_node_point(
      aes(fill = axis, shape = KeyGene_factor, size = is_focal),
      color = "black", stroke = 0.45, alpha = 0.88
    ) +
    geom_node_text(
      aes(label = ifelse(is_focal, name, "")),
      size = 4.6, fontface = "bold", vjust = -1.35
    ) +
    scale_fill_manual(values = pal_axis, na.value = "grey90", name = "Module") +
    scale_shape_manual(
      values = c("TRUE" = 24, "FALSE" = 21),
      name = "Gene type",
      labels = c("KeyGene", "Other")
    ) +
    scale_size_manual(values = c("FALSE" = 3, "TRUE" = 6), guide = "none") +
    scale_edge_width(range = c(0.25, 1.1), guide = "none") +
    theme_graph(base_family = "sans") +
    labs(
      title = paste0(gene, " - 1-hop network: module"),
      subtitle = paste0("correlation threshold >= ", threshold,
                        " | nodes = ", igraph::vcount(graph),
                        " | edges = ", igraph::ecount(graph))
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10.5, hjust = 0.5),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9)
    )

  ggsave(outfile, p, width = 8, height = 7, dpi = 300, bg = "white")
}

plot_deg_network <- function(layout_df, graph, comp, outfile, gene, threshold) {
  labels <- c(
    T25_H1_T25_D1 = "(25C) 1h vs. 1d",
    T42_H1_T42_D1 = "(42C) 1h vs. 1d",
    T25_H1_T42_H1 = "(1h) 25C vs. 42C",
    T25_D1_T42_D1 = "(1d) 25C vs. 42C"
  )

  plot_df <- layout_df
  plot_df$deg_factor <- factor(plot_df[[comp]], levels = c(1, -1))

  p <- ggraph(plot_df) +
    geom_edge_link(aes(width = abs_cor), alpha = 0.28, color = "grey45") +
    geom_node_point(
      aes(fill = deg_factor, shape = KeyGene_factor, size = is_focal),
      color = "black", stroke = 0.45, alpha = 0.9
    ) +
    geom_node_text(
      aes(label = ifelse(is_focal, name, "")),
      size = 4.6, fontface = "bold", vjust = -1.35
    ) +
    scale_fill_manual(
      values = c("1" = "#D73027", "-1" = "#4575B4"),
      labels = c("1" = "Up", "-1" = "Down"),
      na.value = "white",
      name = "DEG"
    ) +
    scale_shape_manual(
      values = c("TRUE" = 24, "FALSE" = 21),
      name = "Gene type",
      labels = c("KeyGene", "Other")
    ) +
    scale_size_manual(values = c("FALSE" = 3, "TRUE" = 6), guide = "none") +
    scale_edge_width(range = c(0.25, 1.1), guide = "none") +
    theme_graph(base_family = "sans") +
    labs(
      title = paste0(gene, " - 1-hop network: DEG pattern"),
      subtitle = paste0(ifelse(comp %in% names(labels), labels[[comp]], comp),
                        " | correlation threshold >= ", threshold,
                        " | nodes = ", igraph::vcount(graph))
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11.5, hjust = 0.5),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9)
    )

  ggsave(outfile, p, width = 8, height = 7, dpi = 300, bg = "white")
}

make_gene_network_outputs <- function(gene, cor_dt, axis_map, node_table, deg_info, output_dir) {
  gene_dir <- file.path(output_dir, "Selected_gene_networks", gene)
  dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

  threshold_info <- choose_threshold(cor_dt, gene, COR_THRESHOLDS, MAX_NEIGHBORS)
  threshold <- threshold_info$threshold[[1]]

  focal_edges <- get_gene_edges_at_threshold(cor_dt, gene, threshold)
  if (!nrow(focal_edges)) {
    warning("No edges found for ", gene, " at threshold ", threshold)
    return(tibble(
      gene = gene,
      chosen_threshold = threshold,
      n_nodes = 0L,
      n_edges = 0L,
      n_focal_neighbors = 0L,
      tried_thresholds = threshold_info$tried_thresholds[[1]],
      output_dir = gene_dir
    ))
  }

  node_ids <- unique(c(gene, focal_edges$to))
  induced_edges <- cor_dt[
    Gene1 %chin% node_ids & Gene2 %chin% node_ids & abs_cor >= threshold,
    .(from = Gene1, to = Gene2, Correlation, abs_cor)
  ]

  g <- graph_from_data_frame(induced_edges, directed = FALSE, vertices = tibble(name = node_ids))

  vertex_df <- tibble(name = V(g)$name) %>%
    left_join(axis_map %>% transmute(name = node, axis, subcluster_id), by = "name") %>%
    left_join(
      node_table %>% select(name, KeyGene, degree, betweenness, closeness, TopCentrality, Bridge),
      by = "name"
    ) %>%
    left_join(deg_info$direction_wide, by = c("name" = "gene")) %>%
    mutate(
      is_focal = name == gene,
      KeyGene = if_else(is.na(KeyGene), FALSE, KeyGene),
      KeyGene_factor = factor(KeyGene, levels = c(TRUE, FALSE)),
      is_DEG_union = name %in% deg_info$deg_union
    )

  for (comp in deg_info$plot_comps) {
    if (!comp %in% colnames(vertex_df)) vertex_df[[comp]] <- 0
    vertex_df[[comp]][is.na(vertex_df[[comp]])] <- 0
  }

  vertex_df <- vertex_df[match(V(g)$name, vertex_df$name), ]
  for (col in colnames(vertex_df)) {
    g <- set_vertex_attr(g, col, value = vertex_df[[col]])
  }

  set.seed(123)
  layout_df <- create_layout(g, layout = "fr")

  edge_detail <- focal_edges %>%
    as_tibble() %>%
    transmute(
      keygene = gene,
      neighbor_gene = to,
      Correlation,
      abs_cor
    ) %>%
    left_join(axis_map %>% transmute(neighbor_gene = node, neighbor_module = axis), by = "neighbor_gene") %>%
    left_join(
      node_table %>% transmute(neighbor_gene = name, neighbor_KeyGene = KeyGene),
      by = "neighbor_gene"
    ) %>%
    mutate(neighbor_DEG_union = neighbor_gene %in% deg_info$deg_union) %>%
    arrange(desc(abs_cor), neighbor_module, neighbor_gene)

  node_info <- vertex_df %>%
    arrange(desc(is_focal), axis, name)

  deg_summary <- bind_rows(lapply(deg_info$plot_comps, function(comp) {
    node_info %>%
      count(comparison = comp, DEG_value = .data[[comp]], name = "n") %>%
      mutate(DEG_status = case_when(
        DEG_value == 1 ~ "Up",
        DEG_value == -1 ~ "Down",
        TRUE ~ "non-DEG"
      ))
  }))

  module_summary <- edge_detail %>%
    count(neighbor_module, name = "n_neighbor") %>%
    arrange(desc(n_neighbor))

  write_csv(sanitize_df(node_info), file.path(gene_dir, paste0(gene, "_1hop_nodes.csv")))
  write_csv(sanitize_df(edge_detail), file.path(gene_dir, paste0(gene, "_1hop_focal_edges.csv")))
  write_csv(sanitize_df(module_summary), file.path(gene_dir, paste0(gene, "_1hop_neighbor_module_summary.csv")))
  write_csv(sanitize_df(deg_summary), file.path(gene_dir, paste0(gene, "_1hop_DEG_pattern_summary.csv")))

  plot_module_network(
    layout_df, g,
    file.path(gene_dir, paste0(gene, "_a_module.png")),
    gene, threshold
  )

  for (comp in deg_info$plot_comps) {
    plot_deg_network(
      layout_df, g, comp,
      file.path(gene_dir, paste0(gene, "_", gsub("_", "-", comp), ".png")),
      gene, threshold
    )
  }

  tibble(
    gene = gene,
    chosen_threshold = threshold,
    n_nodes = vcount(g),
    n_edges = ecount(g),
    n_focal_neighbors = nrow(edge_detail),
    tried_thresholds = threshold_info$tried_thresholds[[1]],
    output_dir = gene_dir
  )
}

## ------------------------------------------------------------------
## 2. DEG union and full-network module inputs
## ------------------------------------------------------------------

deg_info <- read_deg_union_from_rdata(paths$deg_rdata)
write_csv(sanitize_df(deg_info$deg_long), file.path(out_dir, "DEG_union_source_two_heat_comparisons_detail.csv"))
write_csv(tibble(gene = sort(deg_info$deg_union)), file.path(out_dir, "DEG_union_genes_used_for_DEGonly_network.csv"))

message("DEG target comparisons: ", paste(deg_info$target_comps, collapse = ", "))
message("DEG union genes: ", length(deg_info$deg_union))

axis_map <- make_axis_map(paths$node_membership)
write_csv(sanitize_df(axis_map), file.path(out_dir, "full_network_HEAT_ROS_PS_axis_map.csv"))

node_table <- read_csv(paths$full_keygene, show_col_types = FALSE) %>%
  mutate(
    KeyGene = as.logical(KeyGene),
    TopCentrality = as.logical(TopCentrality),
    Bridge = as.logical(Bridge)
  ) %>%
  filter(axis %in% c("HEAT", "ROS", "PS"))

nondeg_keygenes <- node_table %>%
  filter(KeyGene, !(name %in% deg_info$deg_union)) %>%
  transmute(
    keygene = name,
    keygene_module = axis,
    keygene_degree = degree,
    keygene_betweenness = betweenness,
    keygene_closeness = closeness,
    TopCentrality,
    Bridge,
    is_DEG_union = FALSE
  )

write_csv(sanitize_df(nondeg_keygenes), file.path(out_dir, "full_network_nonDEG_KeyGenes_DEGunion_excluded.csv"))
message("Full-network non-DEG KeyGenes in HEAT/ROS/PS: ", nrow(nondeg_keygenes))

## ------------------------------------------------------------------
## 3. Load correlation pairs and build candidate edge table
## ------------------------------------------------------------------

message("Loading full-network correlation pairs...")
cor_dt <- as.data.table(readRDS(paths$cor_pairs))
setnames(cor_dt, old = c("Gene1", "Gene2", "Correlation"), new = c("Gene1", "Gene2", "Correlation"), skip_absent = TRUE)
cor_dt[, abs_cor := abs(Correlation)]

min_threshold <- min(COR_THRESHOLDS)
axis_genes <- unique(axis_map$node)

message("Filtering correlation pairs to HEAT/ROS/PS genes at abs(correlation) >= ", min_threshold, "...")
cor_dt <- cor_dt[
  abs_cor >= min_threshold &
    Gene1 %chin% axis_genes &
    Gene2 %chin% axis_genes,
  .(Gene1, Gene2, Correlation, abs_cor)
]

message("Module-filtered correlation edges: ", nrow(cor_dt))

candidate_genes <- nondeg_keygenes$keygene
edges1 <- cor_dt[Gene1 %chin% candidate_genes,
                 .(keygene = Gene1, neighbor_gene = Gene2, Correlation, abs_cor)]
edges2 <- cor_dt[Gene2 %chin% candidate_genes,
                 .(keygene = Gene2, neighbor_gene = Gene1, Correlation, abs_cor)]
candidate_edges <- rbindlist(list(edges1, edges2), use.names = TRUE) %>%
  as_tibble() %>%
  left_join(nondeg_keygenes %>% select(keygene, keygene_module), by = "keygene") %>%
  left_join(axis_map %>% transmute(neighbor_gene = node, neighbor_module = axis), by = "neighbor_gene") %>%
  filter(!is.na(keygene_module), !is.na(neighbor_module)) %>%
  mutate(
    relation_to_keygene_module = if_else(neighbor_module == keygene_module, "same_module", "other_module"),
    neighbor_is_DEG_union = neighbor_gene %in% deg_info$deg_union
  )

## ------------------------------------------------------------------
## 4. Candidate summary: other two modules and correlation ranking
## ------------------------------------------------------------------

candidate_edges_unique <- candidate_edges %>%
  distinct(keygene, keygene_module, neighbor_gene, neighbor_module, .keep_all = TRUE)

count_wide <- candidate_edges_unique %>%
  count(keygene, keygene_module, neighbor_module, name = "n_neighbor") %>%
  pivot_wider(names_from = neighbor_module, values_from = n_neighbor, values_fill = 0) %>%
  ensure_module_cols()

other_stats <- candidate_edges_unique %>%
  filter(neighbor_module != keygene_module) %>%
  group_by(keygene, keygene_module) %>%
  summarise(
    n_other_module_neighbors = n_distinct(neighbor_gene),
    mean_corr_other = mean(Correlation, na.rm = TRUE),
    mean_abs_corr_other = mean(abs_cor, na.rm = TRUE),
    n_other_DEG_union_neighbors = n_distinct(neighbor_gene[neighbor_is_DEG_union]),
    connected_other_HEAT = collapse_genes(neighbor_gene[neighbor_module == "HEAT"]),
    connected_other_ROS = collapse_genes(neighbor_gene[neighbor_module == "ROS"]),
    connected_other_PS = collapse_genes(neighbor_gene[neighbor_module == "PS"]),
    .groups = "drop"
  )

same_stats <- candidate_edges_unique %>%
  filter(neighbor_module == keygene_module) %>%
  group_by(keygene, keygene_module) %>%
  summarise(
    n_same_module_neighbors = n_distinct(neighbor_gene),
    mean_abs_corr_same = mean(abs_cor, na.rm = TRUE),
    .groups = "drop"
  )

candidate_summary <- count_wide %>%
  left_join(other_stats, by = c("keygene", "keygene_module")) %>%
  left_join(same_stats, by = c("keygene", "keygene_module")) %>%
  left_join(nondeg_keygenes, by = c("keygene", "keygene_module")) %>%
  mutate(
    other_HEAT_count = if_else(keygene_module == "HEAT", NA_integer_, as.integer(HEAT)),
    other_ROS_count = if_else(keygene_module == "ROS", NA_integer_, as.integer(ROS)),
    other_PS_count = if_else(keygene_module == "PS", NA_integer_, as.integer(PS)),
    min_other_module_neighbors = pmin(other_HEAT_count, other_ROS_count, other_PS_count, na.rm = TRUE),
    max_other_module_neighbors = pmax(other_HEAT_count, other_ROS_count, other_PS_count, na.rm = TRUE),
    balance_ratio_other_modules = safe_div(min_other_module_neighbors, max_other_module_neighbors),
    connects_both_other_modules = min_other_module_neighbors > 0
  ) %>%
  mutate(
    across(
      c(n_other_module_neighbors, n_other_DEG_union_neighbors, n_same_module_neighbors),
      ~ replace_na(.x, 0)
    )
  ) %>%
  arrange(
    desc(connects_both_other_modules),
    desc(mean_abs_corr_other),
    desc(balance_ratio_other_modules),
    desc(n_other_module_neighbors),
    keygene
  )

balanced_summary <- candidate_summary %>%
  filter(connects_both_other_modules)

write_csv(sanitize_df(candidate_summary), file.path(out_dir, "nonDEG_KeyGene_cross_module_candidate_summary_all.csv"))
write_csv(sanitize_df(balanced_summary), file.path(out_dir, "nonDEG_KeyGene_connecting_both_other_modules_ranked.csv"))
write_csv(sanitize_df(candidate_edges_unique), file.path(out_dir, "nonDEG_KeyGene_neighbor_edges_at_base_threshold_detail.csv"))

message("Candidates connecting both other modules: ", nrow(balanced_summary))

## ------------------------------------------------------------------
## 5. Choose target genes and draw slide-40-style networks
## ------------------------------------------------------------------

cmd_args <- commandArgs(trailingOnly = TRUE)
genes_arg <- cmd_args[grepl("^--genes=", cmd_args)]
if (length(genes_arg)) {
  TARGET_GENES <- strsplit(sub("^--genes=", "", genes_arg[[1]]), ",", fixed = TRUE)[[1]]
  TARGET_GENES <- trimws(TARGET_GENES)
  TARGET_GENES <- TARGET_GENES[nzchar(TARGET_GENES)]
}

if (!length(TARGET_GENES)) {
  if (isTRUE(PROCESS_ALL_NONDEG_KEYGENES)) {
    TARGET_GENES <- nondeg_keygenes$keygene
  } else {
    TARGET_GENES <- head(balanced_summary$keygene, N_AUTO_GENES)
  }
}

target_tbl <- tibble(gene = TARGET_GENES) %>%
  left_join(candidate_summary, by = c("gene" = "keygene"))

write_csv(sanitize_df(target_tbl), file.path(out_dir, "selected_target_genes.csv"))
message("Selected target genes: ", length(TARGET_GENES))

network_summaries <- bind_rows(lapply(seq_along(TARGET_GENES), function(i) {
  gene <- TARGET_GENES[[i]]
  if (i %% 25 == 0 || i == 1 || i == length(TARGET_GENES)) {
    message("Plotting 1-hop networks: ", i, "/", length(TARGET_GENES), " (", gene, ")")
  }
  tryCatch(
    make_gene_network_outputs(gene, cor_dt, axis_map, node_table, deg_info, out_dir),
    error = function(e) {
      warning("Failed to process ", gene, ": ", conditionMessage(e))
      tibble(
        gene = gene,
        chosen_threshold = NA_real_,
        n_nodes = NA_integer_,
        n_edges = NA_integer_,
        n_focal_neighbors = NA_integer_,
        tried_thresholds = NA_character_,
        output_dir = NA_character_,
        error = conditionMessage(e)
      )
    }
  )
}))

if (nrow(network_summaries)) {
  write_csv(sanitize_df(network_summaries), file.path(out_dir, "selected_target_gene_network_output_summary.csv"))
}

## ------------------------------------------------------------------
## 6. Optional TPM export for selected genes and their 1-hop neighbors
## ------------------------------------------------------------------

selected_node_files <- file.path(out_dir, "Selected_gene_networks", TARGET_GENES, paste0(TARGET_GENES, "_1hop_nodes.csv"))
selected_nodes <- bind_rows(lapply(selected_node_files[file.exists(selected_node_files)], read_csv, show_col_types = FALSE)) %>%
  distinct(name) %>%
  pull(name)

if (file.exists(paths$tpm) && length(selected_nodes)) {
  tpm <- read.table(paths$tpm, header = TRUE, check.names = FALSE, sep = "\t")
  if ("Geneid" %in% colnames(tpm)) {
    tpm_selected <- tpm %>% filter(Geneid %in% selected_nodes)
    write.table(
      tpm_selected,
      file = file.path(out_dir, "selected_target_genes_and_1hop_neighbors_TPM.txt"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
}

## ------------------------------------------------------------------
## 7. Workbook
## ------------------------------------------------------------------

wb <- createWorkbook()
sheets <- list(
  DEG_union = tibble(gene = sort(deg_info$deg_union)),
  nonDEG_KeyGenes = nondeg_keygenes,
  ranked_balanced = balanced_summary,
  all_candidates = candidate_summary,
  selected_targets = target_tbl,
  network_outputs = network_summaries
)

for (nm in names(sheets)) {
  addWorksheet(wb, substr(nm, 1, 31))
  writeData(wb, substr(nm, 1, 31), sanitize_df(sheets[[nm]]))
}

for (sheet in names(wb)) {
  freezePane(wb, sheet, firstRow = TRUE)
  setColWidths(wb, sheet, cols = 1:80, widths = "auto")
}

saveWorkbook(
  wb,
  file.path(out_dir, "Full_nonDEG_KeyGene_balanced_modules_summary.xlsx"),
  overwrite = TRUE
)

cat("\n===== DONE =====\n")
cat("DEG definition: union of ", paste(deg_info$target_comps, collapse = " + "), "\n", sep = "")
cat("DEG union genes:", length(deg_info$deg_union), "\n")
cat("Full-network non-DEG KeyGenes:", nrow(nondeg_keygenes), "\n")
cat("Balanced cross-module candidates:", nrow(balanced_summary), "\n")
cat("Selected target genes:", length(TARGET_GENES), "\n")
cat("Output directory:", out_dir, "\n")
