###############################################################
# 02. Module-level GO-BP ORA
###############################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(clusterProfiler)
  library(org.At.tair.db)
})

set.seed(777)

# === EDIT THESE PATHS FOR YOUR ENVIRONMENT ===
dir.base   <- "/path/to/results"
dir.out    <- file.path(dir.base, "004_Network/module_subnetwork")
file.blast <- "/path/to/Blastp_TAIR11/blastp_araport11_filtered_best.tsv"
# =============================================
if (!dir.exists(dir.out)) dir.create(dir.out, recursive = TRUE)

# Module definitions
source("01_module_definition.R")   # axis_tbl, assign_module()

# Load membership
df_node <- read_csv(
  file.path(dir.base, "004_Network/Leiden_Subgrouping_iterative/<DATE>/Iterative_node_membership_thr200.csv"),
  show_col_types = FALSE
)

# BLASTP best-hit -> TAIR
blast_best <- read_tsv(file.blast, show_col_types = FALSE) %>%
  mutate(geneID = str_replace(qseqid, "\\.t[0-9]+$", ""),
         sseqid = str_replace(sseqid, "\\.[0-9]+$", "")) %>%
  group_by(geneID) %>%
  slice_max(bitscore, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(geneID, sseqid)

# Assign modules and attach TAIR IDs
df_axis <- assign_module(df_node) %>%
  filter(!is.na(axis)) %>%
  left_join(blast_best, by = c("node" = "geneID")) %>%
  dplyr::rename(tair = sseqid)

cat("\n===== Module gene counts =====\n")
print(df_axis %>% group_by(axis) %>%
  summarise(total = n_distinct(node), tair = n_distinct(tair[!is.na(tair)])))

# ── GO-BP ORA ─────────────────────────────────────────────────
axes       <- c("PS", "HEAT", "ROS")
top10_list <- list()

for (ax in axes) {
  cat(sprintf("\n[ %s ] enrichGO BP ...\n", ax))
  genes <- df_axis %>% filter(axis == ax, !is.na(tair)) %>% pull(tair) %>% unique()
  cat("  TAIR:", length(genes), "\n")

  tryCatch({
    er <- enrichGO(gene = genes, OrgDb = org.At.tair.db, keyType = "TAIR",
                   ont = "BP", pAdjustMethod = "BH",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05)

    if (!is.null(er) && nrow(er@result) > 0) {
      res <- as.data.frame(er@result) %>% filter(p.adjust < 0.05) %>% arrange(p.adjust)
      cat("  significant terms:", nrow(res), "| top:", res$Description[1], "\n")

      top10_list[[ax]] <- head(res, 10) %>%
        mutate(axis = ax) %>%
        dplyr::select(axis, ID, Description, p.adjust, Count, GeneRatio)
    } else {
      cat("  no significant terms\n")
    }
  }, error = function(e) cat("  ERROR:", e$message, "\n"))
}

bp_top10 <- bind_rows(top10_list)

# Save
out_csv <- file.path(dir.out, "ORA_PS_HEAT_ROS_BP_top10.csv")
write_csv(bp_top10, out_csv)
cat("\n[Done] saved:", out_csv, "\n")
