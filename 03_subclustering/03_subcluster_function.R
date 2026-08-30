############################################################
## 003. Subcluster ORA (GO Over-Representation Analysis)
############################################################

library(dplyr)
library(readr)
library(stringr)
library(clusterProfiler)
library(org.At.tair.db)
library(tibble)

## ─────────────────────────────────────────────
## 0. Settings
## ─────────────────────────────────────────────
# === PATHS ===
dir.base   <- "/path/to/results"
dir.out    <- "/path/to/results/004_Network/Cluster_function"
file.blast <- "/path/to/Blastp_TAIR11/blastp_araport11_filtered_best.tsv"
file.map   <- file.path(dir.base,
                        "004_Network/Leiden_Subgrouping_iterative/<DATE>/Iterative_node_membership_thr200.csv")
# =============================================

dir.create(dir.out, recursive = TRUE, showWarnings = FALSE)

MIN_TAIR <- 5   # minimum TAIR-mapped genes per subcluster for ORA

## ─────────────────────────────────────────────
## 1. Load membership and BLASTP best-hits
## ─────────────────────────────────────────────
map_df <- read_csv(file.map)

blast_df <- read_tsv(file.blast) %>%
  mutate(
    geneID = str_replace(qseqid, "\\.t[0-9]+$", ""),
    sseqid = str_replace(sseqid, "\\.[0-9]+$",  "")
  )

blast_best <- blast_df %>%
  group_by(geneID) %>%
  slice_max(bitscore, n = 1, with_ties = FALSE) %>%
  ungroup()

# Join TAIR ID onto each Nb gene
map_df <- map_df %>%
  left_join(blast_best %>% dplyr::select(geneID, sseqid),
            by = c("node" = "geneID")) %>%
  dplyr::rename(tair = sseqid)

cat("===== Total subclusters:", length(unique(map_df$subcluster_id)), "=====\n")

## ─────────────────────────────────────────────
## 2. Per-subcluster TAIR counts; select valid subclusters
## ─────────────────────────────────────────────
subcluster_tair_count <- map_df %>%
  filter(!is.na(tair)) %>%
  group_by(subcluster_id) %>%
  summarise(
    total_genes = n(),
    tair_genes  = n_distinct(tair),
    .groups = "drop"
  )

print(head(subcluster_tair_count, 20))

valid_subclusters <- subcluster_tair_count %>%
  filter(tair_genes >= MIN_TAIR) %>%
  pull(subcluster_id)

cat("\n===== Subclusters with >=", MIN_TAIR, "TAIR genes:",
    length(valid_subclusters), "=====\n")

## ─────────────────────────────────────────────
## 3. ORA per subcluster
## ─────────────────────────────────────────────
ora_results_list <- list()
bp_summary_list  <- list()

for (i in seq_along(valid_subclusters)) {
  sub_id <- valid_subclusters[i]
  cat("\n[", i, "/", length(valid_subclusters), "] Subcluster:", sub_id, "\n")

  genes_tair <- map_df %>%
    filter(subcluster_id == !!sub_id, !is.na(tair)) %>%
    pull(tair) %>%
    unique()

  cat("  -> Unique TAIR genes:", length(genes_tair), "\n")

  tryCatch({
    ora_result <- enrichGO(
      gene          = genes_tair,
      OrgDb         = org.At.tair.db,
      keyType       = "TAIR",
      ont           = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05
    )

    if (!is.null(ora_result) && nrow(ora_result@result) > 0) {
      result_df <- as.data.frame(ora_result@result)
      ora_results_list[[sub_id]] <- result_df
      cat("  -> Enriched terms:", nrow(result_df), "\n")

      bp_result <- result_df %>%
        filter(ONTOLOGY == "BP") %>%
        arrange(p.adjust) %>%
        head(1)

      if (nrow(bp_result) > 0) {
        bp_summary_list[[sub_id]] <- data.frame(
          subcluster_id = sub_id,
          ID            = bp_result$ID,
          Description   = bp_result$Description,
          p.adjust      = bp_result$p.adjust,
          Count         = bp_result$Count,
          GeneRatio     = bp_result$GeneRatio,
          stringsAsFactors = FALSE
        )
        cat("  -> Top BP:", bp_result$Description,
            "(p.adj =", format(bp_result$p.adjust, scientific = TRUE, digits = 3), ")\n")
      } else {
        cat("  -> No BP terms\n")
      }
    } else {
      cat("  -> No significant enrichment\n")
      ora_results_list[[sub_id]] <- NULL
    }
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    ora_results_list[[sub_id]] <- NULL
  })
}

## ─────────────────────────────────────────────
## 4. Save outputs
## ─────────────────────────────────────────────
saveRDS(ora_results_list, file.path(dir.out, "ORA_subcluster_results.rds"))

bp_summary_df <- bind_rows(bp_summary_list)
write_csv(bp_summary_df,         file.path(dir.out, "ORA_subcluster_BP_summary.csv"))
write_csv(subcluster_tair_count, file.path(dir.out, "subcluster_tair_info.csv"))

## ─────────────────────────────────────────────
## 5. Summary
## ─────────────────────────────────────────────
cat("\n===== Analysis Summary =====\n")
cat("Total subclusters:", length(unique(map_df$subcluster_id)), "\n")
cat("Valid subclusters (>= ", MIN_TAIR, " TAIR genes):", length(valid_subclusters), "\n", sep = "")
cat("Subclusters with enrichment:",
    length(ora_results_list[!sapply(ora_results_list, is.null)]), "\n")
cat("Subclusters with BP terms:", nrow(bp_summary_df), "\n")

if (nrow(bp_summary_df) > 0) {
  cat("\n===== Top BP Terms by Subcluster (Preview) =====\n")
  print(head(bp_summary_df %>% arrange(p.adjust), 10))
}
