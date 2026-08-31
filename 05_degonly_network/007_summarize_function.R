## Summarize cluster-level functional analysis for new DEG-only network.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(openxlsx)
})

# === PATHS ===
BASE     <- "/path/to/project"
NET      <- file.path(BASE, "Network_DEGonly")
FUNC_DIR <- file.path(NET, "30_Function_Analysis/cluster")
OUT_DIR  <- file.path(BASE, "NetworkComparison_Full_vs_DEGonly")
# ===========================================
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cluster_stat_file          <- file.path(NET, "10_Cluster/Cluster_Stats/<DATE>/Cluster_stats_Leiden_res1.csv")
subcluster_membership_file <- file.path(NET, "20_SubCluster/Leiden_Subgrouping_iterative/<DATE>/Iterative_node_membership_thr200.csv")
all_ora_file <- file.path(FUNC_DIR, "ORA_cluster_ALL_results.csv")

out_csv <- file.path(OUT_DIR, "DEGonly_cluster_function_top_terms.csv")
out_xlsx <- file.path(OUT_DIR, "DEGonly_cluster_function_top_terms.xlsx")

cluster_stats <- read_csv(cluster_stat_file, show_col_types = FALSE) %>%
  mutate(cluster_id = paste0("C", cluster)) %>%
  select(cluster, cluster_id, node_count, edge_count)

subcluster_counts <- read_csv(subcluster_membership_file, show_col_types = FALSE) %>%
  mutate(cluster = as.integer(res1_cluster)) %>%
  summarise(subcluster_count = n_distinct(subcluster_id), .by = cluster)

ora <- read_csv(all_ora_file, show_col_types = FALSE) %>%
  mutate(
    cluster_num = as.integer(str_remove(cluster, "^C")),
    minus_log10_padj = -log10(p.adjust)
  )

top_by_ontology_long <- ora %>%
  group_by(cluster_num, ONTOLOGY) %>%
  arrange(p.adjust, pvalue, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    cluster = cluster_num,
    ONTOLOGY,
    top_ID = ID,
    top_Description = Description,
    top_GeneRatio = GeneRatio,
    top_Count = Count,
    top_p.adjust = p.adjust,
    top_minus_log10_padj = minus_log10_padj
  )

top_by_ontology_wide <- top_by_ontology_long %>%
  mutate(ONTOLOGY = recode(ONTOLOGY, BP = "BP", CC = "CC", MF = "MF")) %>%
  pivot_wider(
    names_from = ONTOLOGY,
    values_from = c(top_ID, top_Description, top_GeneRatio, top_Count, top_p.adjust, top_minus_log10_padj),
    names_glue = "{ONTOLOGY}_{.value}"
  )

top_overall <- ora %>%
  group_by(cluster_num) %>%
  arrange(p.adjust, pvalue, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    cluster = cluster_num,
    top_overall_ontology = ONTOLOGY,
    top_overall_ID = ID,
    top_overall_Description = Description,
    top_overall_GeneRatio = GeneRatio,
    top_overall_Count = Count,
    top_overall_p.adjust = p.adjust,
    top_overall_minus_log10_padj = minus_log10_padj
  )

summary <- cluster_stats %>%
  left_join(subcluster_counts, by = "cluster") %>%
  left_join(top_overall, by = "cluster") %>%
  left_join(top_by_ontology_wide, by = "cluster") %>%
  arrange(cluster)

table_for_paper <- summary %>%
  transmute(
    cluster,
    node_count,
    edge_count,
    subcluster_count,
    top_enriched_ontology = top_overall_ontology,
    top_enriched_GO = top_overall_ID,
    top_enriched_term = top_overall_Description,
    top_enriched_p.adjust = top_overall_p.adjust,
    BP_top_term = BP_top_Description,
    BP_p.adjust = BP_top_p.adjust,
    CC_top_term = CC_top_Description,
    CC_p.adjust = CC_top_p.adjust,
    MF_top_term = MF_top_Description,
    MF_p.adjust = MF_top_p.adjust
  )

top5_by_ontology <- ora %>%
  group_by(cluster_num, ONTOLOGY) %>%
  arrange(p.adjust, pvalue, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  transmute(
    cluster = cluster_num,
    ONTOLOGY,
    ID,
    Description,
    GeneRatio,
    Count,
    p.adjust,
    minus_log10_padj
  ) %>%
  arrange(cluster, ONTOLOGY, p.adjust)

write_csv(table_for_paper, out_csv)

wb <- createWorkbook()
add_sheet <- function(sheet, data) {
  data <- data %>%
    mutate(across(where(is.character), ~ iconv(.x, from = "", to = "UTF-8", sub = " ")))
  addWorksheet(wb, sheet)
  writeData(wb, sheet, data)
  freezePane(wb, sheet, firstRow = TRUE)
  addFilter(wb, sheet, rows = 1, cols = seq_len(ncol(data)))
  setColWidths(wb, sheet, cols = seq_len(ncol(data)), widths = "auto")
}

add_sheet("Table_for_paper", table_for_paper)
add_sheet("Top_overall", top_overall)
add_sheet("Top_by_ontology", top_by_ontology_long)
add_sheet("Top5_by_ontology", top5_by_ontology)
add_sheet("All_ORA_results", ora)
saveWorkbook(wb, out_xlsx, overwrite = TRUE)

message("Wrote: ", out_csv)
message("Wrote: ", out_xlsx)
