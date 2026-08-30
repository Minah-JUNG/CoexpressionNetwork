##------------------------------------------------
## 2025-10-15  Leiden vs Louvain 
##------------------------------------------------

rm(list = ls())

library(dplyr)
library(readr)
library(igraph)
library(tidygraph)
library(stringr)
library(lubridate)
library(tibble)

##------------------------------------------------
## Settings
##------------------------------------------------
set.seed(12345)   # reproducibility

# === PATH ===
base_dir <- "/path/to/results/004_Network"
# ===========================================
setwd(base_dir)

today <- Sys.Date()

make_dir <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

out_dir_base     <- file.path("./Network_bases", today)
out_dir_leiden   <- file.path("./Leiden_Results", today)
out_dir_compare  <- file.path("./Compare", today)
out_dir_node     <- file.path("./Table_NodeInfo", today)
out_dir_cluster  <- file.path("./Cluster_Stats", today)

sapply(list(out_dir_base, out_dir_leiden, out_dir_compare, 
            out_dir_node, out_dir_cluster), make_dir)

log_file <- file.path(out_dir_base, "Leiden_vs_Louvain_runlog.txt")
cat("Run log\n", file = log_file)

script_start <- Sys.time()
cat("Script start:", script_start, "\n", file = log_file, append = TRUE)

##------------------------------------------------
## Load correlation data
##------------------------------------------------
file_corr <- "../003_Correlation/Filtered_correlation_pairs.rds"

DataCorr <- readRDS(file_corr) %>%
  filter(abs(Correlation) >= 0.7, Pvalue < 0.05) %>%
  rename(Gene_1 = 1, Gene_2 = 2)

##------------------------------------------------
## Remove edges between two singleton nodes
##------------------------------------------------
deg_temp <- table(c(DataCorr$Gene_1, DataCorr$Gene_2))
single_nodes <- names(deg_temp[deg_temp == 1])

DataCorr_f <- DataCorr %>%
  filter(!(Gene_1 %in% single_nodes & Gene_2 %in% single_nodes))

cat("Original pairs:", nrow(DataCorr),
    "Filtered pairs:", nrow(DataCorr_f), "\n",
    file = log_file, append = TRUE)

##------------------------------------------------
## Build graph
##------------------------------------------------
g <- graph_from_data_frame(DataCorr_f[, c("Gene_1", "Gene_2")], directed = FALSE)
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
g <- delete_vertices(g, which(degree(g) == 0))

cat("Graph stats - Node:", vcount(g), 
    "Edge:", ecount(g), "\n",
    file = log_file, append = TRUE)

##------------------------------------------------
## Compute centralities
##------------------------------------------------
deg_vec <- degree(g)
btw_vec <- betweenness(g, normalized = TRUE)
cls_vec <- closeness(g, normalized = TRUE)

##------------------------------------------------
## Leiden sweep
##------------------------------------------------
res_values <- seq(1.0, 3.0, 0.25)

leiden_results <- list()
leiden_summary <- list()

for (res in res_values) {
  
  cat("\n----- Leiden (resolution =", res, ") -----\n", 
      file = log_file, append = TRUE)
  
  tt <- system.time({
    cl <- cluster_leiden(g, resolution = res, objective_function = "modularity")
  })
  
  # 저장
  saveRDS(cl, file.path(out_dir_leiden, paste0("Leiden_res", res, "_comm.rds")))
  
  memb <- membership(cl)
  memb_df <- tibble(Gene = names(memb), cluster = as.integer(memb))
  write_csv(memb_df, file.path(out_dir_node, paste0("Leiden_res", res, "_membership.csv")))
  
  cluster_stats <- tibble(
    cluster = names(sizes(cl)),
    size = as.integer(sizes(cl))
  ) %>% arrange(desc(size))
  
  write_csv(cluster_stats, file.path(out_dir_cluster, paste0("Leiden_res", res, "_cluster_stats.csv")))
  
  # tbl_graph 저장
  g_tg_tmp <- as_tbl_graph(g) %>%
    activate(nodes) %>%
    mutate(
      cluster = memb[name],
      Degree = deg_vec[name],
      Betweenness = btw_vec[name],
      Closeness = cls_vec[name]
    )
  saveRDS(g_tg_tmp, file.path(out_dir_leiden, paste0("Global_g_tg_Leiden_res", res, ".rds")))
  
  # summary row 저장
  leiden_results[[as.character(res)]] <- list(comm = cl)
  leiden_summary[[as.character(res)]] <- tibble(
    resolution = res,
    time_sec = tt["elapsed"],
    n_clusters = length(cl),
    top5 = paste(head(sort(sizes(cl), TRUE), 5), collapse = ", ")
  )
}

leiden_summary_df <- bind_rows(leiden_summary)
write_csv(leiden_summary_df, file.path(out_dir_leiden, "Leiden_summary_modularity.csv"))

##------------------------------------------------
## Louvain
##------------------------------------------------
tt_lv <- system.time({
  cl_louvain <- cluster_louvain(g)
})
saveRDS(cl_louvain, file.path(out_dir_base, "Global_louvain.rds"))

##------------------------------------------------
## Compare Leiden vs Louvain
##------------------------------------------------
compare_list <- lapply(res_values, function(res) {
  cl_leiden <- leiden_results[[as.character(res)]]$comm
  ari <- tryCatch(compare(cl_louvain, cl_leiden, method = "adjusted.rand"), error = function(e) NA)
  nmi <- tryCatch(compare(cl_louvain, cl_leiden, method = "nmi"), error = function(e) NA)
  
  tibble(
    resolution = res,
    louvain_clusters = length(cl_louvain),
    leiden_clusters = length(cl_leiden),
    ari = ari,
    nmi = nmi
  )
})

compare_df <- bind_rows(compare_list)
write_csv(compare_df, file.path(out_dir_compare, "Leiden_vs_Louvain_modularity_compare.csv"))

##------------------------------------------------
## Membership table (Gene | Louvain | Leiden_* )
##------------------------------------------------
membership_df <- tibble(Gene = V(g)$name,
                        louvain = as.integer(membership(cl_louvain)[V(g)$name]))

for (res in res_values) {
  cl <- leiden_results[[as.character(res)]]$comm
  membership_df[[paste0("leiden_", formatC(res, format = "f", digits = 2))]] <-
    as.integer(membership(cl)[membership_df$Gene])
}

write_csv(membership_df, file.path(out_dir_node, "Gene_membership_table_leiden_louvain.csv"))

##------------------------------------------------
## Auto-select best Leiden resolution (target = 10 clusters)
##------------------------------------------------
target_k <- length(cl_louvain)
leiden_summary_df$abs_diff <- abs(leiden_summary_df$n_clusters - target_k)

chosen <- leiden_summary_df %>%
  arrange(abs_diff, time_sec) %>% slice(1)

chosen_res <- chosen$resolution
chosen_comm <- leiden_results[[as.character(chosen_res)]]$comm

##------------------------------------------------
## Find bridge genes
##------------------------------------------------
find_bridge_genes <- function(graph, memb) {
  ed <- as.data.frame(as_edgelist(graph), FALSE)
  colnames(ed) <- c("from", "to")
  ed <- ed %>%
    mutate(from_cl = memb[from], to_cl = memb[to]) %>%
    filter(from_cl != to_cl)
  unique(c(ed$from, ed$to))
}

bridge_genes <- find_bridge_genes(g, membership(chosen_comm))

##------------------------------------------------
## Node annotation (hub/connector/bridge/key)
##------------------------------------------------
node_tbl <- tibble(
  name = V(g)$name,
  cluster_leiden = membership(chosen_comm)[V(g)$name],
  degree = deg_vec[V(g)$name],
  betweenness = btw_vec[V(g)$name],
  closeness = cls_vec[V(g)$name]
)

thr_deg <- quantile(node_tbl$degree, 0.95)
thr_btw <- quantile(node_tbl$betweenness, 0.95)
thr_cls <- quantile(node_tbl$closeness, 0.95)

node_tbl <- node_tbl %>%
  mutate(
    HubGene = ifelse(degree >= thr_deg, "HubGene", "Other"),
    ConnectorGene = ifelse(betweenness >= thr_btw | closeness >= thr_cls, "ConnectorGene", "Other"),
    BridgeGene = ifelse(name %in% bridge_genes, "BridgeGene", "Other"),
    KeyGene = ifelse(
      (HubGene == "HubGene" | ConnectorGene == "ConnectorGene") &
        BridgeGene == "BridgeGene",
      "KeyGene", "Other")
  )

write_csv(node_tbl, file.path(out_dir_node, paste0("Node_annotation_chosenLeiden_res", chosen_res, ".csv")))

##------------------------------------------------
## Cluster stats (chosen Leiden)
##------------------------------------------------
cluster_stats2 <- lapply(sort(unique(node_tbl$cluster_leiden)), function(cid) {
  vids <- V(g)[membership(chosen_comm) == cid]
  sg <- induced_subgraph(g, vids)
  tibble(cluster = cid, node_count = vcount(sg), edge_count = ecount(sg))
})

cluster_stats2_df <- bind_rows(cluster_stats2)
write_csv(cluster_stats2_df, file.path(out_dir_cluster, paste0("Cluster_Stats_ChosenLeiden_res", chosen_res, ".csv")))

##------------------------------------------------
## Save tbl_graph (chosen Leiden)
##------------------------------------------------
g_tg <- as_tbl_graph(g) %>%
  activate(nodes) %>%
  mutate(
    cluster = membership(chosen_comm)[name],
    Degree = deg_vec[name],
    Betweenness = btw_vec[name],
    Closeness = cls_vec[name],
    HubGene = node_tbl$HubGene[match(name, node_tbl$name)],
    ConnectorGene = node_tbl$ConnectorGene[match(name, node_tbl$name)],
    BridgeGene = node_tbl$BridgeGene[match(name, node_tbl$name)],
    KeyGene = node_tbl$KeyGene[match(name, node_tbl$name)]
  )

saveRDS(g_tg, file.path(out_dir_leiden, paste0("Global_g_tg_Leiden_chosen_res", chosen_res, ".rds")))

##------------------------------------------------
## Log 종료
##------------------------------------------------
script_end <- Sys.time()

cat("Script end:", script_end, "\n", file = log_file, append = TRUE)
cat("Total elapsed (mins):", 
    difftime(script_end, script_start, units = "mins"), "\n",
    file = log_file, append = TRUE)
cat("Auto-chosen resolution:", chosen_res,
    "clusters:", chosen$n_clusters, "\n",
    file = log_file, append = TRUE)

cat("All done.\nOutput dirs:\n",
    out_dir_leiden, "\n", out_dir_base, "\n", out_dir_node, "\n",
    out_dir_cluster, "\n", out_dir_compare, "\n",
    "Log:", log_file, "\n")
