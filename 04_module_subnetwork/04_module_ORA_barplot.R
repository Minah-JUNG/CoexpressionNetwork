###############################################################
# GO BP Top-10 barplot — PS / HEAT / ROS modules
###############################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(clusterProfiler)
  library(org.At.tair.db)
})

set.seed(777)

# === PATHS ===
dir.base   <- "/path/to/results"
dir.out    <- file.path(dir.base, "004_Network/module_subnetwork")
file.blast <- "/path/to/Blastp_TAIR11/blastp_araport11_filtered_best.tsv"
# =============================================
if (!dir.exists(dir.out)) dir.create(dir.out, recursive = TRUE)

# Module definitions
source("01_module_definition.R")   # axis_tbl, axis_palette, assign_module()

# Load membership
df_node <- read_csv(
  file.path(dir.base,
            "004_Network/Leiden_Subgrouping_iterative/<DATE>/Iterative_node_membership_thr200.csv"),
  show_col_types = FALSE
)

# BLASTP best-hit -> TAIR
blast_best <- read_tsv(file.blast, show_col_types = FALSE) %>%
  mutate(
    geneID = str_replace(qseqid, "\\.t[0-9]+$", ""),
    sseqid = str_replace(sseqid, "\\.[0-9]+$", "")
  ) %>%
  group_by(geneID) %>%
  slice_max(bitscore, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(geneID, sseqid)

# Assign modules + attach TAIR
df_axis <- assign_module(df_node) %>%
  filter(!is.na(axis)) %>%
  left_join(blast_best, by = c("node" = "geneID")) %>%
  rename(tair = sseqid)

cat("\n===== Module gene counts =====\n")
print(
  df_axis %>%
    group_by(axis) %>%
    summarise(
      total_Nb  = n_distinct(node),
      tair_hits = n_distinct(tair[!is.na(tair)])
    )
)

# ── 4. GO-BP ORA ────────────────────────────────────────────────
axes        <- c("PS", "HEAT", "ROS")
ora_results <- list()
top10_list  <- list()

for (ax in axes) {
  cat(sprintf("\n[ %-4s ] enrichGO BP ...\n", ax))

  gene_vec <- df_axis %>%
    filter(axis == ax, !is.na(tair)) %>%
    pull(tair) %>% unique()

  cat("  TAIR genes:", length(gene_vec), "\n")

  tryCatch({
    er <- enrichGO(
      gene          = gene_vec,
      OrgDb         = org.At.tair.db,
      keyType       = "TAIR",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = FALSE
    )

    if (!is.null(er) && nrow(er@result) > 0) {
      res <- as.data.frame(er@result) %>%
        filter(p.adjust < 0.05) %>%
        arrange(p.adjust)

      ora_results[[ax]] <- res
      cat("  Significant BP terms:", nrow(res), "\n")
      cat("  Top term:", res$Description[1], "\n")

      top10_list[[ax]] <- head(res, 10) %>%
        mutate(axis = ax) %>%
        dplyr::select(axis, ID, Description, p.adjust, Count, GeneRatio)
    } else {
      cat("  No significant enrichment\n")
    }
  }, error = function(e) cat("  ERROR:", e$message, "\n"))
}

bp_top10 <- bind_rows(top10_list)

# ── 5. Save ORA results ──────────────────────────────────────────
saveRDS(ora_results, file.path(dir.out, "ORA_PS_HEAT_ROS_BP_results.rds"))
write_csv(bp_top10,  file.path(dir.out, "ORA_PS_HEAT_ROS_BP_top10.csv"))
cat("\n[saved] ORA_PS_HEAT_ROS_BP_results.rds\n")
cat("[saved] ORA_PS_HEAT_ROS_BP_top10.csv\n")

# ── 6. Color palette (from Multiple_Cluster_network.png) ─────────
#   PS   : medium green   (#5CB85C)
#   HEAT : coral-red pink (#E8736C)
#   ROS  : golden yellow  (#F5C842)
module_colors <- c(
  PS   = "#5CB85C",
  HEAT = "#E8736C",
  ROS  = "#F5C842"
)
module_colors_bg <- sapply(module_colors, function(x) adjustcolor(x, alpha.f = 0.18))

module_labels <- c(
  PS   = "PS — Energy Production & Maintenance",
  HEAT = "HEAT — Proteostasis & Cellular Protection",
  ROS  = "ROS — Signal Integration & Defense"
)

# ── 7. Plot settings (function_top5 스타일) ──────────────────────
base_size       <- 16
desc_size       <- 4.7
id_size         <- 4.0
title_size      <- 22
axis_title_size <- 14

# ── 8. Prepare plot data ────────────────────────────────────────
df_plot <- bp_top10 %>%
  mutate(
    axis          = factor(axis, levels = axes),
    log10p        = -log10(p.adjust),
    GeneRatio_num = sapply(GeneRatio, function(x) eval(parse(text = x)))
  ) %>%
  group_by(axis) %>%
  arrange(p.adjust) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# 세 module 공통 x축 최대값
xmax_val <- max(df_plot$log10p, na.rm = TRUE) * 1.15

# ── 9. Single panel function ────────────────────────────────────
make_panel <- function(df_ax, axis_code, xmax) {

  # y축: p.adjust 오름차순 → 위에서 아래로 (ggplot reversed)
  df_ax <- df_ax %>% arrange(p.adjust) %>%
    mutate(term_id = factor(ID, levels = rev(ID)))

  col_fg <- module_colors[axis_code]
  col_bg <- module_colors_bg[axis_code]

  ggplot(df_ax, aes(x = log10p, y = term_id)) +
    # ── 배경 막대 (full-width) ──────────────────────────────
    geom_col(aes(x = xmax), fill = col_bg, width = 0.82) +
    # ── 실제 값 막대 ────────────────────────────────────────
    geom_col(fill = col_fg, color = NA, width = 0.82) +
    # ── GO Description 텍스트 (막대 안, 왼쪽 정렬) ──────────
    geom_text(
      aes(x = 0.025 * xmax, label = Description),
      hjust  = 0, vjust = 0.25,
      size   = desc_size,
      color  = "grey10",
      fontface = "plain"
    ) +
    # ── GO ID 텍스트 (Description 위에 살짝) ────────────────
    geom_text(
      aes(x = 0.025 * xmax, label = ID),
      hjust  = 0, vjust = 2.0,
      size   = id_size,
      color  = "grey40"
    ) +
    scale_x_continuous(
      limits = c(0, xmax),
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(
      title = module_labels[axis_code],
      x     = expression(-log[10]~(p.adjust)),
      y     = NULL
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      plot.title         = element_text(face = "bold", size = 15,
                                        hjust = 0, margin = margin(b = 6)),
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      axis.title.x       = element_text(size = axis_title_size),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
      plot.margin        = margin(10, 18, 10, 15)
    )
}

# ── 10. Build panels ────────────────────────────────────────────
panel_list <- lapply(axes, function(ax) {
  make_panel(
    df_plot %>% filter(axis == ax) %>% arrange(p.adjust),
    ax, xmax_val
  )
})
names(panel_list) <- axes

# ── 11. Legend ──────────────────────────────────────────────────
legend_df <- data.frame(
  module = factor(axes, levels = axes), x = 1, y = 1
)
p_leg <- ggplot(legend_df, aes(x = x, y = y, color = module)) +
  geom_point(size = 5) +
  scale_color_manual(
    values = module_colors,
    name   = "Module",
    labels = module_labels,
    guide  = guide_legend(nrow = 1, title.position = "left")
  ) +
  theme_void() +
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_text(face = "bold", size = 17),
    legend.text      = element_text(size = 13),
    legend.key       = element_rect(fill = "white", color = NA),
    legend.key.size  = unit(1.1, "lines")
  )
legend_grob <- cowplot::get_legend(p_leg)

# ── 12. Combine & save (가로 3패널) ─────────────────────────────
final <- wrap_plots(panel_list, ncol = 3) +   # ← 가로 정렬
  plot_annotation(
    title = "Gene Ontology — Top 10 BP Terms per Module",
    theme = theme(
      plot.title = element_text(face = "bold", size = title_size,
                                hjust = 0.5, margin = margin(b = 10))
    )
  )

final_with_legend <- cowplot::plot_grid(
  final, legend_grob,
  ncol        = 1,
  rel_heights = c(1, 0.06)   # legend 비율 살짝 늘림 (가로 레이아웃)
)

out_png <- file.path(dir.out, "GO_BP_top10_PS_HEAT_ROS.png")
out_pdf <- file.path(dir.out, "GO_BP_top10_PS_HEAT_ROS.pdf")

ggsave(out_png, final_with_legend, width = 18, height = 8, dpi = 300, bg = "white")  # 가로
ggsave(out_pdf, final_with_legend, width = 18, height = 8,             bg = "white")

cat("\n[saved]", out_png, "\n")
cat("[saved]", out_pdf, "\n")
cat("\n===== ALL DONE =====\n")
