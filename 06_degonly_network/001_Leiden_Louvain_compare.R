##================================================
## Gene Co-expression Network — Leiden vs Louvain 비교
## Fixed seed: 12345
##================================================

rm(list = ls())

SEED <- 12345
set.seed(SEED)

library(dplyr)
library(readr)
library(igraph)
library(tidygraph)
library(tibble)

## ─────────────────────────────────────────────
## 0. Settings
## ─────────────────────────────────────────────
FILE_CORR   <- "../Filtered_correlation_pairs.rds"
CORR_CUTOFF <- 0.7
P_CUTOFF    <- 0.05
RES_VALUES  <- seq(1.0, 3.0, by = 0.25)

today <- format(Sys.Date(), "%Y-%m-%d")
make_dir <- function(path) { dir.create(path, recursive = TRUE, showWarnings = FALSE); path }

out_leiden  <- make_dir(paste0("./Leiden_Results/",  today, "/"))
out_compare <- make_dir(paste0("./Compare/",          today, "/"))
out_cluster <- make_dir(paste0("./Cluster_Stats/",    today, "/"))
log_file    <- paste0(out_leiden, "LeidenLouvain_compare.log")

log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
cat("", file = log_file, append = FALSE)
log_msg("=== Script start: ", Sys.time(), " ===")
log_msg("seed fixed: ", SEED)

## ─────────────────────────────────────────────
## 1. Load & filter
## ─────────────────────────────────────────────
DataCorr <- readRDS(FILE_CORR) %>%
  setNames(c("Gene_1", "Gene_2", "Correlation", "Pvalue")) %>%
  filter(abs(Correlation) >= CORR_CUTOFF, Pvalue < P_CUTOFF)
log_msg("Loaded pairs: ", nrow(DataCorr))

deg_raw    <- table(c(DataCorr$Gene_1, DataCorr$Gene_2))
singletons <- names(deg_raw[deg_raw == 1])
DataCorr   <- DataCorr %>%
  filter(!(Gene_1 %in% singletons & Gene_2 %in% singletons))
log_msg("After singleton removal: ", nrow(DataCorr))

## ─────────────────────────────────────────────
## 2. Build graph
## ─────────────────────────────────────────────
g <- graph_from_data_frame(DataCorr[, c("Gene_1", "Gene_2")], directed = FALSE) %>%
  simplify(remove.multiple = TRUE, remove.loops = TRUE) %>%
  delete_vertices(which(degree(.) == 0))

log_msg("Graph — nodes: ", vcount(g), " | edges: ", ecount(g))

## ─────────────────────────────────────────────
## 3. Louvain (1회)
## ─────────────────────────────────────────────
log_msg("Running Louvain...")
set.seed(SEED)
tt_lv <- system.time(cl_louvain <- cluster_louvain(g))
log_msg("Louvain done — clusters: ", length(cl_louvain),
        " | time: ", round(tt_lv["elapsed"], 1), "s")

saveRDS(cl_louvain, paste0(out_leiden, "Louvain_comm.rds"))

write_csv(
  tibble(cluster = names(sizes(cl_louvain)), size = as.integer(sizes(cl_louvain))) %>%
    arrange(desc(size)),
  paste0(out_cluster, "Louvain_cluster_stats.csv")
)

## ─────────────────────────────────────────────
## 4. Leiden sweep
## ─────────────────────────────────────────────
leiden_results <- list()
leiden_summary <- list()

for (resv in RES_VALUES) {
  log_msg("----- Leiden resolution = ", resv, " -----")

  set.seed(SEED)
  tt <- system.time(
    cl <- cluster_leiden(g, resolution = resv, objective_function = "modularity")
  )
  log_msg("  clusters: ", length(cl), " | time: ", round(tt["elapsed"], 1), "s")

  saveRDS(cl, paste0(out_leiden, "Leiden_res", resv, "_comm.rds"))

  write_csv(
    tibble(cluster = names(sizes(cl)), size = as.integer(sizes(cl))) %>%
      arrange(desc(size)),
    paste0(out_cluster, "Leiden_res", resv, "_cluster_stats.csv")
  )

  leiden_results[[as.character(resv)]] <- list(comm = cl, time = tt["elapsed"])
  leiden_summary[[as.character(resv)]] <- tibble(
    resolution = resv,
    n_clusters = length(cl),
    top5       = paste(head(sort(sizes(cl), decreasing = TRUE), 5), collapse = ", "),
    time_sec   = round(tt["elapsed"], 1)
  )
}

summary_df <- bind_rows(leiden_summary)
write_csv(summary_df, paste0(out_leiden, "Leiden_summary_modularity.csv"))

## ─────────────────────────────────────────────
## 5. Leiden vs Louvain 비교
## ─────────────────────────────────────────────
compare_df <- lapply(RES_VALUES, function(resv) {
  cl <- leiden_results[[as.character(resv)]]$comm
  ari <- tryCatch(compare(cl_louvain, cl, method = "adjusted.rand"), error = function(e) NA)
  nmi <- tryCatch(compare(cl_louvain, cl, method = "nmi"),           error = function(e) NA)
  tibble(
    resolution       = resv,
    louvain_clusters = length(cl_louvain),
    leiden_clusters  = length(cl),
    ARI              = round(ari, 3),
    NMI              = round(nmi, 3),
    distance_to_1    = round(sqrt((1 - ari)^2 + (1 - nmi)^2), 3),
    score            = round((ari + nmi) / 2, 3)
  )
}) %>% bind_rows()

write_csv(compare_df, paste0(out_compare, "Leiden_vs_Louvain_compare.csv"))

chosen_df <- compare_df %>%
  filter(!is.na(distance_to_1)) %>%
  arrange(distance_to_1, desc(score), desc(ARI), desc(NMI), resolution) %>%
  slice(1)

write_csv(chosen_df, paste0(out_compare, "Chosen_resolution_by_ARI_NMI.csv"))
log_msg("Chosen resolution closest to ARI/NMI = 1: ", chosen_df$resolution,
        " | Leiden clusters: ", chosen_df$leiden_clusters,
        " | ARI: ", chosen_df$ARI,
        " | NMI: ", chosen_df$NMI,
        " | distance_to_1: ", chosen_df$distance_to_1,
        " | score: ", chosen_df$score)

log_msg("=== Script end: ", Sys.time(), " ===")
cat("\nComparison result:\n")
print(compare_df)
cat("\nChosen resolution closest to ARI/NMI = 1:\n")
print(chosen_df)
cat("\nDone. Check:", out_compare, "\n")
