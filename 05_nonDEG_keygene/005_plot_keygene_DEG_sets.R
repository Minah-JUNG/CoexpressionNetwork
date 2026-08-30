###############################################################
# Full subnetwork KeyGene DEG-set comparison
###############################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
  library(tibble)
})

# === PATHS ===
out_root <- "/path/to/output_root"
out_dir  <- file.path(out_root,
  paste0("Full_subnetwork_keygene_DEG_set_comparison_", format(Sys.Date(), "%Y-%m-%d")))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

paths <- list(
  full_keygene_nodes = "/path/to/Full_vs_DEGonly_subnetwork/full_subnetwork_node_centrality_KeyGene.csv",
  deg_binary         = "/path/to/DEG/DEG_matrix_binary.txt",
  deg_rdata          = "/path/to/DEG/Nicotiana_DEGs_Heat.RData"
)
# =============================================

comparisons_set1 <- c("T25_H1_T42_H1", "T25_D1_T42_D1")
comparisons_set2 <- c("T25_H1_T42_H1", "T25_D1_T42_D1", "T42_H1_T42_D1", "T25_H1_T25_D1")

comparison_labels <- c(
  T25_H1_T42_H1 = "1h: 25C vs 42C",
  T25_D1_T42_D1 = "1d: 25C vs 42C",
  T42_H1_T42_D1 = "42C: 1h vs 1d",
  T25_H1_T25_D1 = "25C: 1h vs 1d"
)

read_deg_binary <- function(path) {
  comparison_cols <- strsplit(readLines(path, n = 1, warn = FALSE), "\t", fixed = TRUE)[[1]]
  read_tsv(
    path,
    skip = 1,
    col_names = c("gene", comparison_cols),
    show_col_types = FALSE
  ) %>%
    mutate(across(-gene, ~ suppressWarnings(as.numeric(.x))))
}

load_deseq_results <- function(path) {
  env <- new.env()
  load(path, envir = env)
  if (!exists("result_list", envir = env)) {
    stop("result_list not found in: ", path)
  }

  bind_rows(lapply(names(env$result_list), function(comp) {
    as.data.frame(env$result_list[[comp]]) %>%
      rownames_to_column("gene") %>%
      as_tibble() %>%
      transmute(
        gene,
        comparison = comp,
        baseMean,
        log2FoldChange,
        foldChange = 2^log2FoldChange,
        abs_log2FoldChange = abs(log2FoldChange),
        abs_foldChange = 2^abs(log2FoldChange),
        lfcSE,
        stat,
        pvalue,
        padj
      )
  }))
}

summarise_de_stats <- function(de_res, comparisons, suffix) {
  de_res %>%
    filter(comparison %in% comparisons) %>%
    group_by(gene) %>%
    summarise(
      "max_abs_log2FC_{suffix}" := max(abs_log2FoldChange, na.rm = TRUE),
      "max_abs_FC_{suffix}" := max(abs_foldChange, na.rm = TRUE),
      "min_pvalue_{suffix}" := suppressWarnings(min(pvalue, na.rm = TRUE)),
      "min_padj_{suffix}" := suppressWarnings(min(padj, na.rm = TRUE)),
      "comparison_max_abs_log2FC_{suffix}" := comparison[which.max(abs_log2FoldChange)][1],
      .groups = "drop"
    ) %>%
    mutate(
      across(starts_with("max_abs"), ~ if_else(is.infinite(.x), NA_real_, .x)),
      across(starts_with("min_"), ~ if_else(is.infinite(.x), NA_real_, .x))
    )
}

plot_pie <- function(df, set_name, outfile) {
  pie_df <- df %>%
    count(.data[[paste0("DEG_", set_name)]], name = "n") %>%
    mutate(
      DEG_status = if_else(.data[[paste0("DEG_", set_name)]], "DEG", "non-DEG"),
      pct = n / sum(n) * 100,
      label = paste0(DEG_status, "\n", n, " (", sprintf("%.1f", pct), "%)")
    )

  p <- ggplot(pie_df, aes(x = "", y = n, fill = DEG_status)) +
    geom_col(width = 1, color = "white", linewidth = 0.5) +
    coord_polar(theta = "y") +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4) +
    scale_fill_manual(values = c("DEG" = "#D95F02", "non-DEG" = "#4D4D4D")) +
    labs(title = paste0("Full Subnetwork KeyGenes: ", set_name), fill = NULL) +
    theme_void(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom"
    )

  ggsave(outfile, p, width = 5.5, height = 5.2, dpi = 300)
  pie_df
}

plot_distribution <- function(non_deg_long, value_col, x_label, title, outfile, log_x = FALSE) {
  p <- ggplot(non_deg_long, aes(x = .data[[value_col]], fill = DEG_set)) +
    geom_density(alpha = 0.35, color = NA, adjust = 1.1) +
    geom_rug(alpha = 0.18, sides = "b") +
    facet_wrap(~ DEG_set, ncol = 1, scales = "free_y") +
    labs(title = title, x = x_label, y = "Density") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none",
      strip.background = element_rect(fill = "grey92", color = "grey75"),
      strip.text = element_text(face = "bold")
    )

  if (log_x) {
    p <- p + scale_x_log10()
  }

  ggsave(outfile, p, width = 7.2, height = 6.2, dpi = 300)
  p
}

keygenes <- read_csv(paths$full_keygene_nodes, show_col_types = FALSE) %>%
  filter(KeyGene) %>%
  select(name, module = axis, degree, betweenness, closeness, TopCentrality, Bridge, KeyGene)

deg_binary <- read_deg_binary(paths$deg_binary)
de_res <- load_deseq_results(paths$deg_rdata)

deg_flags <- deg_binary %>%
  transmute(
    name = gene,
    DEG_set1 = rowSums(across(all_of(comparisons_set1)), na.rm = TRUE) > 0,
    DEG_set2 = rowSums(across(all_of(comparisons_set2)), na.rm = TRUE) > 0
  )

stats_set1 <- summarise_de_stats(de_res, comparisons_set1, "set1") %>%
  rename(name = gene)
stats_set2 <- summarise_de_stats(de_res, comparisons_set2, "set2") %>%
  rename(name = gene)

comparison_detail <- keygenes %>%
  left_join(deg_flags, by = "name") %>%
  left_join(stats_set1, by = "name") %>%
  left_join(stats_set2, by = "name") %>%
  mutate(
    DEG_set1 = if_else(is.na(DEG_set1), FALSE, DEG_set1),
    DEG_set2 = if_else(is.na(DEG_set2), FALSE, DEG_set2),
    DEG_class = case_when(
      DEG_set1 ~ "DEG_set1",
      DEG_set2 & !DEG_set1 ~ "DEG_set2_only",
      TRUE ~ "non_DEG_set2"
    )
  ) %>%
  arrange(DEG_class, module, desc(degree), name)

summary_tbl <- bind_rows(
  comparison_detail %>%
    count(DEG_set1, name = "n_keygenes") %>%
    mutate(DEG_set = "set1_first_two_union", DEG_status = if_else(DEG_set1, "DEG", "non-DEG")) %>%
    select(DEG_set, DEG_status, n_keygenes),
  comparison_detail %>%
    count(DEG_set2, name = "n_keygenes") %>%
    mutate(DEG_set = "set2_all_four_union", DEG_status = if_else(DEG_set2, "DEG", "non-DEG")) %>%
    select(DEG_set, DEG_status, n_keygenes)
) %>%
  group_by(DEG_set) %>%
  mutate(percent = n_keygenes / sum(n_keygenes) * 100) %>%
  ungroup()

summary_by_module <- bind_rows(
  comparison_detail %>%
    count(module, DEG_set1, name = "n_keygenes") %>%
    mutate(DEG_set = "set1_first_two_union", DEG_status = if_else(DEG_set1, "DEG", "non-DEG")) %>%
    select(DEG_set, module, DEG_status, n_keygenes),
  comparison_detail %>%
    count(module, DEG_set2, name = "n_keygenes") %>%
    mutate(DEG_set = "set2_all_four_union", DEG_status = if_else(DEG_set2, "DEG", "non-DEG")) %>%
    select(DEG_set, module, DEG_status, n_keygenes)
) %>%
  group_by(DEG_set, module) %>%
  mutate(percent = n_keygenes / sum(n_keygenes) * 100) %>%
  ungroup()

non_deg_set1 <- comparison_detail %>%
  filter(!DEG_set1) %>%
  transmute(
    name, module,
    DEG_set = "set1_non_DEG",
    max_abs_log2FC = max_abs_log2FC_set1,
    max_abs_FC = max_abs_FC_set1,
    min_pvalue = min_pvalue_set1,
    min_padj = min_padj_set1,
    neg_log10_min_pvalue = -log10(min_pvalue),
    neg_log10_min_padj = -log10(min_padj)
  )

non_deg_set2 <- comparison_detail %>%
  filter(!DEG_set2) %>%
  transmute(
    name, module,
    DEG_set = "set2_non_DEG",
    max_abs_log2FC = max_abs_log2FC_set2,
    max_abs_FC = max_abs_FC_set2,
    min_pvalue = min_pvalue_set2,
    min_padj = min_padj_set2,
    neg_log10_min_pvalue = -log10(min_pvalue),
    neg_log10_min_padj = -log10(min_padj)
  )

non_deg_distribution <- bind_rows(non_deg_set1, non_deg_set2)

dist_summary <- non_deg_distribution %>%
  group_by(DEG_set) %>%
  summarise(
    n_non_DEG_keygenes = n(),
    median_max_abs_log2FC = median(max_abs_log2FC, na.rm = TRUE),
    q75_max_abs_log2FC = quantile(max_abs_log2FC, 0.75, na.rm = TRUE),
    max_max_abs_log2FC = max(max_abs_log2FC, na.rm = TRUE),
    median_min_pvalue = median(min_pvalue, na.rm = TRUE),
    q25_min_pvalue = quantile(min_pvalue, 0.25, na.rm = TRUE),
    median_min_padj = median(min_padj, na.rm = TRUE),
    q25_min_padj = quantile(min_padj, 0.25, na.rm = TRUE),
    .groups = "drop"
  )

pie_set1 <- plot_pie(
  comparison_detail,
  "set1",
  file.path(out_dir, "pie_full_keygenes_DEG_set1_first_two_union.png")
)
pie_set2 <- plot_pie(
  comparison_detail,
  "set2",
  file.path(out_dir, "pie_full_keygenes_DEG_set2_all_four_union.png")
)

plot_distribution(
  non_deg_distribution,
  "max_abs_log2FC",
  "max |log2 fold change| among set comparisons",
  "Non-DEG Full KeyGenes: Fold-change Distribution",
  file.path(out_dir, "density_nonDEG_keygenes_max_abs_log2FC_set1_vs_set2.png")
)

plot_distribution(
  non_deg_distribution %>% filter(is.finite(neg_log10_min_pvalue)),
  "neg_log10_min_pvalue",
  "-log10(min pvalue) among set comparisons",
  "Non-DEG Full KeyGenes: p-value Distribution",
  file.path(out_dir, "density_nonDEG_keygenes_neglog10_min_pvalue_set1_vs_set2.png")
)

plot_distribution(
  non_deg_distribution %>% filter(is.finite(neg_log10_min_padj)),
  "neg_log10_min_padj",
  "-log10(min adjusted pvalue) among set comparisons",
  "Non-DEG Full KeyGenes: adjusted p-value Distribution",
  file.path(out_dir, "density_nonDEG_keygenes_neglog10_min_padj_set1_vs_set2.png")
)

comparison_long <- non_deg_distribution %>%
  select(name, module, DEG_set, max_abs_log2FC, neg_log10_min_pvalue, neg_log10_min_padj) %>%
  pivot_longer(
    cols = c(max_abs_log2FC, neg_log10_min_pvalue, neg_log10_min_padj),
    names_to = "metric",
    values_to = "value"
  )

p_box <- ggplot(comparison_long, aes(x = DEG_set, y = value, fill = DEG_set)) +
  geom_violin(alpha = 0.35, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.18, outlier.alpha = 0.25) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  labs(
    title = "Non-DEG Full KeyGenes: Fold-change and p-value Metrics",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    strip.background = element_rect(fill = "grey92", color = "grey75"),
    strip.text = element_text(face = "bold")
  )
ggsave(file.path(out_dir, "violin_nonDEG_keygenes_metrics_set1_vs_set2.png"), p_box, width = 7, height = 8, dpi = 300)

write_csv(comparison_detail, file.path(out_dir, "full_keygene_DEG_set_comparison_detail.csv"))
write_csv(summary_tbl, file.path(out_dir, "full_keygene_DEG_set_summary.csv"))
write_csv(summary_by_module, file.path(out_dir, "full_keygene_DEG_set_summary_by_module.csv"))
write_csv(non_deg_distribution, file.path(out_dir, "full_keygene_nonDEG_distribution_metrics.csv"))
write_csv(dist_summary, file.path(out_dir, "full_keygene_nonDEG_distribution_summary.csv"))

wb <- createWorkbook()
addWorksheet(wb, "DEG_summary")
writeData(wb, "DEG_summary", summary_tbl)
addWorksheet(wb, "DEG_summary_by_module")
writeData(wb, "DEG_summary_by_module", summary_by_module)
addWorksheet(wb, "nonDEG_dist_summary")
writeData(wb, "nonDEG_dist_summary", dist_summary)
addWorksheet(wb, "keygene_detail")
writeData(wb, "keygene_detail", comparison_detail)
addWorksheet(wb, "nonDEG_distribution")
writeData(wb, "nonDEG_distribution", non_deg_distribution)
for (sheet in names(wb)) {
  freezePane(wb, sheet, firstRow = TRUE)
  setColWidths(wb, sheet, cols = 1:50, widths = "auto")
}
saveWorkbook(wb, file.path(out_dir, "Full_subnetwork_keygene_DEG_set_comparison.xlsx"), overwrite = TRUE)

cat("\n===== DONE =====\n")
cat("Output directory:", out_dir, "\n\n")
cat("[DEG summary]\n")
print(summary_tbl)
cat("\n[Non-DEG distribution summary]\n")
print(dist_summary)
