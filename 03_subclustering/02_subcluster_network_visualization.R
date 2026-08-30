############################################################
## 002. Subcluster Network Visualization
############################################################

rm(list = ls())
set.seed(12345)

library(igraph)
library(dplyr)
library(readr)
library(ggraph)
library(ggplot2)

## ─────────────────────────────────────────────
## 0. Settings
## ─────────────────────────────────────────────
# === PATHS ===
file.map     <- "/path/to/results/004_Network/Leiden_Subgrouping_iterative/<DATE>/Iterative_node_membership_thr200.csv"
file.corr    <- "/path/to/results/003_Correlation/Filtered_correlation_pairs.rds"
dir.out.base <- "/path/to/results/004_Network/SubclusterNetworks"
# =============================================

# Correlation thresholds at which to draw each subcluster network.
th.cor.list <- c(0.7, 0.9)

# Minimum component size to keep in the "filtered" graph.
MIN_COMP    <- 3

# Optional: restrict to a subset of subclusters. Leave NULL to draw all.
# Example: target_subclusters <- c("C_1_4_2_1_1", "C_6_3")
target_subclusters <- NULL

## ─────────────────────────────────────────────
## 1. Load data
## ─────────────────────────────────────────────
map_df  <- read_csv(file.map)
corr_df <- readRDS(file.corr)

sub_list <- unique(map_df$subcluster_id)
if (!is.null(target_subclusters)) {
  sub_list <- intersect(sub_list, target_subclusters)
}
cat("Subclusters to draw:", length(sub_list), "\n")

dir.create(dir.out.base, recursive = TRUE, showWarnings = FALSE)

## ─────────────────────────────────────────────
## 2. Per-threshold, per-subcluster loop
## ─────────────────────────────────────────────
for (th.cor in th.cor.list) {

  message("\n==============================================")
  message("Correlation threshold: ", th.cor)
  message("==============================================\n")

  dir.thr <- file.path(dir.out.base, paste0("thr_", th.cor))
  dir.create(dir.thr, recursive = TRUE, showWarnings = FALSE)

  corr_df_th <- corr_df %>% filter(abs(Correlation) > th.cor)

  for (target_sub in sub_list) {

    message("  Subcluster: ", target_sub)
    dir.sub <- file.path(dir.thr, target_sub)
    dir.create(dir.sub, showWarnings = FALSE)

    # Nodes belonging to this subcluster
    target_nodes <- map_df %>%
      filter(subcluster_id == target_sub) %>%
      pull(node) %>%
      unique()

    write_csv(data.frame(node = target_nodes),
              file.path(dir.sub, "nodes.csv"))

    # Edges within the subcluster
    sub_edges <- corr_df_th %>%
      filter(Gene1 %in% target_nodes & Gene2 %in% target_nodes)

    write_csv(sub_edges, file.path(dir.sub, "edges.csv"))

    ## --- Raw graph
    g_sub <- graph_from_data_frame(
      sub_edges,
      directed = FALSE,
      vertices = data.frame(name = target_nodes)
    )
    saveRDS(g_sub, file.path(dir.sub, "graph_raw.rds"))

    p_raw <- ggraph(g_sub, layout = "graphopt") +
      geom_edge_link(aes(alpha = abs(Correlation)), edge_colour = "grey70") +
      geom_node_point(size = 3, color = "steelblue") +
      geom_node_text(aes(label = name), size = 2.5, repel = TRUE) +
      theme_void() +
      ggtitle(paste0(target_sub, " (raw, n=", length(target_nodes), ")"))

    ggsave(file.path(dir.sub, "Network_raw.png"),
           p_raw, width = 12, height = 10, dpi = 300)

    ## --- Filtered graph: drop degree-0 nodes; keep components of size >= MIN_COMP
    g_filtered <- delete_vertices(g_sub, which(degree(g_sub) == 0))

    comp <- components(g_filtered)
    keep_nodes <- names(comp$membership)[comp$csize[comp$membership] >= MIN_COMP]

    if (length(keep_nodes) == 0) {
      message("     Filtered graph: no component >= ", MIN_COMP, " nodes. Skip.")
      next
    }

    g_final <- induced_subgraph(g_filtered, vids = keep_nodes)
    saveRDS(g_final, file.path(dir.sub, "graph_filtered.rds"))

    p_final <- ggraph(g_final, layout = "fr") +
      geom_edge_link(aes(alpha = abs(Correlation)), edge_colour = "grey70") +
      geom_node_point(size = 3.5, color = "steelblue") +
      geom_node_text(aes(label = name), size = 2.5, repel = TRUE) +
      theme_void() +
      ggtitle(paste0(target_sub, " (filtered, n=", gorder(g_final), ")"))

    ggsave(file.path(dir.sub, "Network_filtered.png"),
           p_final, width = 12, height = 10, dpi = 300)

    message("     Completed.")
  }

  message("\nThreshold ", th.cor, " completed.\n")
}

cat("\nAll thresholds completed. Output under:\n  ", dir.out.base, "\n")
