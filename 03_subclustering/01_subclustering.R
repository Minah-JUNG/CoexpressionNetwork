##================================================
## 001. Iterative Leiden Subclustering
##================================================

rm(list = ls())
set.seed(12345)

library(igraph)
library(tidygraph)
library(dplyr)
library(readr)
library(stringr)

## ─────────────────────────────────────────────
## 0. Settings
## ─────────────────────────────────────────────
# === PATHS ===
FILE_GRAPH <- "/path/to/results/004_Network/Leiden_Results/<DATE>/Global_g_tg_Leiden_chosen_res1.rds"

today   <- format(Sys.Date(), "%Y-%m-%d")
out_dir <- paste0("/path/to/results/004_Network/Leiden_Subgrouping_iterative/", today, "/")
# =============================================
sub_dir <- paste0(out_dir, "Subgraphs_thr200/")
dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)

SIZE_THR  <- 200
MIN_NODES <- 3
RES_BASE  <- 1

## ─────────────────────────────────────────────
## 1. Load global graph
## ─────────────────────────────────────────────
cat("===== Loading graph...\n")
g_tg    <- readRDS(FILE_GRAPH)
g       <- as.igraph(g_tg)
node_df <- g_tg %>% activate(nodes) %>% as_tibble()

cat("===== Graph loaded:", vcount(g), "nodes,", ecount(g), "edges\n")

## ─────────────────────────────────────────────
## 2. Iterative Leiden function
## ─────────────────────────────────────────────
iterative_leiden <- function(graph, parent_name, threshold, res, min_nodes) {
  results       <- list()
  cluster_queue <- list(list(graph = graph, parent = parent_name))

  while (length(cluster_queue) > 0) {
    current       <- cluster_queue[[1]]
    cluster_queue <- cluster_queue[-1]
    g_sub         <- current$graph
    pname         <- current$parent

    if (vcount(g_sub) <= threshold) {
      if (vcount(g_sub) >= min_nodes) {
        results <- c(results, list(list(
          cluster_id = pname,
          nodes      = V(g_sub)$name,
          size       = vcount(g_sub),
          edges      = ecount(g_sub)
        )))
      }
      next
    }

    cl   <- cluster_leiden(g_sub, resolution = res, objective_function = "modularity")
    memb <- membership(cl)

    for (cid in sort(unique(memb))) {
      vids <- V(g_sub)[memb == cid]
      subg <- induced_subgraph(g_sub, vids)
      n    <- vcount(subg)
      if (n < min_nodes) next

      child_name <- paste0(pname, "_", cid)

      if (n > threshold) {
        cluster_queue <- c(cluster_queue,
                           list(list(graph = subg, parent = child_name)))
      } else {
        results <- c(results, list(list(
          cluster_id = child_name,
          nodes      = V(subg)$name,
          size       = n,
          edges      = ecount(subg)
        )))
      }
    }
  }
  return(results)
}

## ─────────────────────────────────────────────
## 3. Run iterative Leiden per global cluster
## ─────────────────────────────────────────────
all_results <- list()

for (cid in sort(unique(node_df$cluster))) {
  nodes_in <- node_df$name[node_df$cluster == cid]
  subg     <- induced_subgraph(g, vids = nodes_in)
  n        <- vcount(subg)
  pname    <- paste0("C_", cid)   # "C_1", "C_2", ... -> downstream-compatible

  cat(">>> Global cluster", cid, "| nodes:", n, "-> ")

  if (n < MIN_NODES) {
    cat("too small, skip\n")
    next
  }

  if (n <= SIZE_THR) {
    cat("saved as-is\n")
    all_results <- c(all_results, list(list(
      cluster_id = pname,
      nodes      = V(subg)$name,
      size       = n,
      edges      = ecount(subg)
    )))
    next
  }

  sub_results <- iterative_leiden(
    graph       = subg,
    parent_name = pname,
    threshold   = SIZE_THR,
    res         = RES_BASE,
    min_nodes   = MIN_NODES
  )
  cat(length(sub_results), "subclusters\n")
  all_results <- c(all_results, sub_results)
}

cat("\n===== Total subclusters:", length(all_results), "\n")

## ─────────────────────────────────────────────
## 4. Save results
## ─────────────────────────────────────────────

## (1) cluster stats
cluster_stats <- lapply(all_results, function(x) {
  data.frame(cluster_id = x$cluster_id,
             node_count = x$size,
             edge_count = x$edges,
             stringsAsFactors = FALSE)
}) %>% bind_rows()

write_csv(cluster_stats,
          paste0(out_dir, "Iterative_cluster_stats_thr", SIZE_THR, ".csv"))

## (2) node membership
node_membership <- lapply(all_results, function(x) {
  parent_part <- str_extract(x$cluster_id, "(?<=^C_)\\d+")
  data.frame(node           = x$nodes,
             parent_cluster = parent_part,
             subcluster_id  = x$cluster_id,
             stringsAsFactors = FALSE)
}) %>% bind_rows()

write_csv(node_membership,
          paste0(out_dir, "Iterative_node_membership_thr", SIZE_THR, ".csv"))

## (3) summary stats
summary_stats <- cluster_stats %>%
  summarise(
    total_clusters = n(),
    mean_size      = mean(node_count),
    median_size    = median(node_count),
    max_size       = max(node_count),
    min_size       = min(node_count)
  )

write_csv(summary_stats,
          paste0(out_dir, "Iterative_summary_thr", SIZE_THR, ".csv"))

## (4) subgraph RDS + full results
cat("===== Saving subgraph RDS files...\n")
for (x in all_results) {
  subg <- induced_subgraph(g, vids = x$nodes)
  saveRDS(subg,
          file = paste0(sub_dir, "Subcluster_", x$cluster_id, "_thr", SIZE_THR, ".rds"))
}
saveRDS(all_results,
        paste0(out_dir, "Iterative_iterresults_thr", SIZE_THR, "_full.rds"))

cat("===== Saved", length(all_results), "subcluster .rds files to", sub_dir, "\n")

## ─────────────────────────────────────────────
## 5. sessionInfo for reproducibility
## ─────────────────────────────────────────────
session_file <- paste0(out_dir, "sessionInfo_", format(Sys.Date(), "%Y%m%d"), ".txt")
sink(session_file)
cat("==============================================\n")
cat("Session Info for Reproducibility\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("seed: 12345\n")
cat("==============================================\n\n")
print(sessionInfo())
sink()

cat("\nAll completed. Results saved under:\n", out_dir, "\n")
