##--------------------------------------------------
## Iterative Leiden Clustering — res1 클러스터 기준 시작
## Fixed seed: 12345
##--------------------------------------------------

rm(list = ls())
set.seed(12345)

library(igraph)
library(tidygraph)
library(dplyr)
library(readr)

## ─────────────────────────────────────────────
## 0. Settings
## ─────────────────────────────────────────────
gnetwork_day <- "<DATE>"
# === PATH ===
base_dir     <- "/path/to/Network_DEGonly/20_SubCluster"
# ===========================================

FILE_GRAPH   <- paste0(base_dir, "/../10_Cluster/Leiden_Results/", gnetwork_day,
                       "/g_tg_Leiden_FINAL_res1.rds")

SIZE_THR  <- 200   # 이 이하가 될 때까지 쪼갬
MIN_NODES <- 3     # 이보다 작은 subcluster는 버림
RES_BASE  <- 1     # leiden resolution

out_dir  <- paste0(base_dir, "/Leiden_Subgrouping_iterative/", gnetwork_day, "/")
sub_dir  <- paste0(out_dir, "Subgraphs_thr", SIZE_THR, "/")
dir.create(sub_dir,  recursive = TRUE, showWarnings = FALSE)

## ─────────────────────────────────────────────
## 1. Load global graph
## ─────────────────────────────────────────────
cat("Loading graph...\n")
g_tbl <- readRDS(FILE_GRAPH)
g     <- as.igraph(g_tbl)   # tbl_graph → igraph
cat("Graph loaded:", vcount(g), "nodes,", ecount(g), "edges\n")

## ─────────────────────────────────────────────
## 2. Extract res1 cluster subgraphs
## ─────────────────────────────────────────────
# tbl_graph의 node 테이블에서 cluster 컬럼 읽기
node_df <- g_tbl %>% activate(nodes) %>% as_tibble()
cat("Columns in node table:", paste(colnames(node_df), collapse = ", "), "\n")

# cluster 컬럼명 확인 (pipeline에서 'cluster'로 저장)
res1_clusters <- sort(unique(node_df$cluster))
cat("res1 clusters found:", length(res1_clusters), "\n\n")

## ─────────────────────────────────────────────
## 3. Iterative Leiden function
## ─────────────────────────────────────────────
iterative_leiden <- function(graph, parent_name, threshold, res, min_nodes) {
  results      <- list()
  cluster_queue <- list(list(graph = graph, parent = parent_name))

  while (length(cluster_queue) > 0) {
    current       <- cluster_queue[[1]]
    cluster_queue <- cluster_queue[-1]
    g_sub         <- current$graph
    pname         <- current$parent

    # 이미 threshold 이하면 바로 저장
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

    # threshold 초과 → leiden으로 쪼개기
    cl   <- cluster_leiden(g_sub, resolution = res, objective_function = "modularity")
    memb <- membership(cl)

    for (cid in sort(unique(memb))) {
      vids  <- V(g_sub)[memb == cid]
      subg  <- induced_subgraph(g_sub, vids)
      n     <- vcount(subg)
      if (n < min_nodes) next

      child_name <- paste0(pname, "_", cid)

      if (n > threshold) {
        # 아직 크면 queue에 추가
        cluster_queue <- c(cluster_queue, list(list(graph = subg, parent = child_name)))
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
## 4. Main: res1 클러스터별로 iterative leiden 실행
## ─────────────────────────────────────────────
all_results <- list()

for (cid in res1_clusters) {
  cat("▶ Processing res1 cluster", cid, "...\n")

  # res1 cluster subgraph 추출
  nodes_in_cluster <- node_df$name[node_df$cluster == cid]
  subg_res1 <- induced_subgraph(g, vids = nodes_in_cluster)
  n <- vcount(subg_res1)
  cat("  nodes:", n, "\n")

  parent_name <- paste0("C", cid)

  # threshold 이하면 그냥 그대로 저장
  if (n <= SIZE_THR) {
    if (n >= MIN_NODES) {
      all_results <- c(all_results, list(list(
        cluster_id = parent_name,
        nodes      = V(subg_res1)$name,
        size       = n,
        edges      = ecount(subg_res1)
      )))
    }
    cat("  → Small cluster, saved as-is.\n")
    next
  }

  # threshold 초과 → iterative leiden
  sub_results <- iterative_leiden(
    graph     = subg_res1,
    parent_name = parent_name,
    threshold = SIZE_THR,
    res       = RES_BASE,
    min_nodes = MIN_NODES
  )
  cat("  → Split into", length(sub_results), "subclusters.\n")
  all_results <- c(all_results, sub_results)
}

cat("\nTotal subclusters:", length(all_results), "\n")

## ─────────────────────────────────────────────
## 5. Save results
## ─────────────────────────────────────────────

## (1) Cluster stats
cluster_stats <- lapply(all_results, function(x) {
  tibble(cluster_id = x$cluster_id, node_count = x$size, edge_count = x$edges)
}) %>% bind_rows() %>% arrange(cluster_id)

write_csv(cluster_stats, paste0(out_dir, "Iterative_cluster_stats_thr", SIZE_THR, ".csv"))

## (2) Node membership
node_membership <- lapply(all_results, function(x) {
  # res1 상위 클러스터 번호 추출 (C1, C2, ...)
  top_cluster <- sub("^C(\\d+).*", "\\1", x$cluster_id)
  tibble(
    node           = x$nodes,
    res1_cluster   = as.integer(top_cluster),
    subcluster_id  = x$cluster_id
  )
}) %>% bind_rows()

write_csv(node_membership, paste0(out_dir, "Iterative_node_membership_thr", SIZE_THR, ".csv"))

## (3) Summary stats
summary_stats <- cluster_stats %>%
  summarise(
    total_clusters = n(),
    mean_size      = round(mean(node_count), 1),
    median_size    = median(node_count),
    max_size       = max(node_count),
    min_size       = min(node_count)
  )
write_csv(summary_stats, paste0(out_dir, "Iterative_summary_thr", SIZE_THR, ".csv"))

## (4) Subgraph RDS files
cat("Saving subgraph RDS files...\n")
for (x in all_results) {
  subg <- induced_subgraph(g, vids = x$nodes)
  saveRDS(subg, paste0(sub_dir, "Subcluster_", x$cluster_id, "_thr", SIZE_THR, ".rds"))
}

## (5) Full results list
saveRDS(all_results, paste0(out_dir, "Iterative_iterresults_thr", SIZE_THR, "_full.rds"))

cat("\nDone. Results saved to:\n", out_dir, "\n")
print(summary_stats)
