###############################################################
# SET1 final panel plots for keygene DEG / non-DEG comparison
# a. KeyGene DEG proportion
# b. DEG-keygene connected genes: DEG status
# c. DEG-keygene connected genes: module with DEG status
# d. DEG-keygene connected DEG genes: function keywords
# e. DEG-keygene connected non-DEG genes: function keywords
# f. non-DEG-keygene connected genes: DEG status
# g. non-DEG-keygene connected genes: module with DEG status
# h. non-DEG-keygene connected DEG genes: function keywords
# i. non-DEG-keygene connected non-DEG genes: function keywords
###############################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(grid)
})

# === PATHS ===
out_root <- "/path/to/output_root"
out_dir  <- file.path(out_root,
  paste0("SET1_KeyGene_DEG_nonDEG_final_panels_", format(Sys.Date(), "%Y-%m-%d")))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  keygene_detail     = "/path/to/SET1_keygene_DEG_detail.csv",
  degKG_neighbors    = "/path/to/SET1_DEG_keygene_unique_connected_genes_annotation.csv",
  nondegKG_neighbors = "/path/to/SET1_unique_connected_genes_function_annotation.csv"
)
# =============================================

status_colors_main <- c("DEG" = "#C43C2B", "non-DEG" = "#8DB9D8")
status_colors_red <- c("DEG" = "#C43C2B", "non-DEG" = "#F4B6A6")
status_colors_blue <- c("DEG" = "#2F6FA8", "non-DEG" = "#B7D4EA")

term_colors <- c(
  rna_translation = "#4776C6",
  heat_chaperone = "#F47C2C",
  signaling_kinase = "#9E9E9E",
  metabolism = "#FFC20A",
  stress_defense = "#6FA8DC",
  hormone = "#5FA85B",
  photosynthesis_chloroplast = "#2F5597",
  ros_redox = "#A85A16",
  transcription = "#5B5B5B",
  transport = "#9A8F21",
  protein_turnover = "#1F4E79"
)

term_labels <- c(
  rna_translation = "RNA\ntranslation",
  heat_chaperone = "heat\nchaperone",
  signaling_kinase = "signaling\nkinase",
  metabolism = "metabolism",
  stress_defense = "stress\ndefense",
  hormone = "hormone",
  photosynthesis_chloroplast = "photosynthesis\nchloroplast",
  ros_redox = "ROS\nredox",
  transcription = "transcription",
  transport = "transport",
  protein_turnover = "protein\nturnover"
)

keyword_patterns <- list(
  photosynthesis_chloroplast = "photosynth|chloroplast|thylakoid|photosystem|light-harvesting|light harvesting|chlorophyll",
  heat_chaperone = "heat shock|chaperone|HSP|DnaJ|heat stress",
  ros_redox = "reactive oxygen|ROS|redox|peroxidase|oxidase|thioredoxin|glutathione|superoxide|catalase",
  stress_defense = "stress|defense|pathogen|immune|disease|wound|hypersensitive",
  hormone = "abscisic|ABA|auxin|ethylene|jasmon|salicyl|gibberellin|cytokinin|brassinosteroid",
  transcription = "transcription|DNA-binding|DNA binding|transcription factor|zinc finger|MYB|bHLH|WRKY|NAC",
  signaling_kinase = "kinase|phosphatase|calcium|calmodulin|receptor|signaling|signal transduction",
  transport = "transport|transporter|channel|carrier|ABC transporter|aquaporin",
  protein_turnover = "ubiquitin|proteasome|E3 ligase|protein degradation|autophagy",
  rna_translation = "ribosomal|ribosome|translation|RNA|splicing|spliceosome|tRNA|rRNA",
  metabolism = "metabolic|metabolism|biosynthesis|catabolic|synthase|dehydrogenase|transferase"
)

count_keywords <- function(df) {
  text <- paste(
    ifelse(is.na(df$short_description), "", df$short_description),
    ifelse(is.na(df$Curator_summary), "", df$Curator_summary),
    ifelse(is.na(df$Computational_description), "", df$Computational_description)
  )

  bind_rows(lapply(names(keyword_patterns), function(category) {
    tibble(
      category = category,
      term = term_labels[[category]],
      n = sum(grepl(keyword_patterns[[category]], text, ignore.case = TRUE), na.rm = TRUE)
    )
  })) %>%
    filter(n > 0) %>%
    arrange(desc(n), category)
}

plot_status_pie <- function(df, outfile, title, colors, show_percent = FALSE) {
  plot_df <- df %>%
    mutate(
      status = factor(status, levels = c("DEG", "non-DEG")),
      percent = n / sum(n) * 100
    )
  if (show_percent) {
    plot_df <- plot_df %>%
      mutate(label = paste0(status, "\n", n, " (", sprintf("%.1f", percent), "%)"))
  } else {
    plot_df <- plot_df %>%
      mutate(label = paste0(status, "\n", n))
  }

  title_element <- if (is.null(title) || identical(title, "")) {
    element_blank()
  } else {
    element_text(face = "bold", hjust = 0.5, size = 16)
  }

  p <- ggplot(plot_df, aes(x = "", y = n, fill = status)) +
    geom_col(width = 1, color = "white", linewidth = 0.8) +
    coord_polar(theta = "y") +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 6, fontface = "bold") +
    scale_fill_manual(values = colors, drop = FALSE) +
    labs(title = title, fill = NULL) +
    theme_void(base_size = 16) +
    theme(
      plot.title = title_element,
      legend.position = "bottom",
      legend.text = element_text(size = 12, face = "bold")
    )

  ggsave(outfile, p, width = 5.8, height = 5.4, dpi = 300, bg = "white")
}

plot_module_status_bar <- function(df, outfile, title, colors) {
  plot_df <- df %>%
    mutate(
      module = factor(module, levels = c("HEAT", "PS", "ROS")),
      status = factor(status, levels = c("DEG", "non-DEG"))
    )

  title_element <- if (is.null(title) || identical(title, "")) {
    element_blank()
  } else {
    element_text(face = "bold", hjust = 0.5, size = 16)
  }

  p <- ggplot(plot_df, aes(x = module, y = n, fill = status)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.5) +
    geom_text(
      aes(label = n),
      position = position_stack(vjust = 0.5),
      size = 5.2,
      fontface = "bold",
      color = "black"
    ) +
    scale_fill_manual(values = colors, drop = FALSE) +
    labs(title = title, x = NULL, y = "Gene count", fill = NULL) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = title_element,
      axis.text = element_text(size = 13, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      legend.position = "right",
      legend.text = element_text(size = 12, face = "bold")
    )

  ggsave(outfile, p, width = 6.4, height = 4.8, dpi = 300, bg = "white")
}

binary_treemap <- function(df, x = 0, y = 0, w = 1, h = 1) {
  if (nrow(df) == 1) {
    return(df %>% mutate(xmin = x, xmax = x + w, ymin = y, ymax = y + h))
  }

  total <- sum(df$n)
  cumulative <- cumsum(df$n)
  split_at <- which.min(abs(cumulative - total / 2))
  split_at <- max(1, min(split_at, nrow(df) - 1))

  left <- df[seq_len(split_at), , drop = FALSE]
  right <- df[(split_at + 1):nrow(df), , drop = FALSE]
  left_fraction <- sum(left$n) / total

  if (w >= h) {
    w_left <- w * left_fraction
    bind_rows(
      binary_treemap(left, x, y, w_left, h),
      binary_treemap(right, x + w_left, y, w - w_left, h)
    )
  } else {
    h_bottom <- h * left_fraction
    bind_rows(
      binary_treemap(left, x, y, w, h_bottom),
      binary_treemap(right, x, y + h_bottom, w, h - h_bottom)
    )
  }
}

plot_keyword_treemap <- function(summary_df, outfile, title) {
  plot_df <- summary_df %>%
    arrange(desc(n), category)
  rect_df <- binary_treemap(plot_df) %>%
    mutate(
      cx = (xmin + xmax) / 2,
      cy = (ymin + ymax) / 2,
      area = (xmax - xmin) * (ymax - ymin),
      label = if_else(area >= 0.035, paste0(term, "\n", n), ""),
      text_size = pmax(3.6, pmin(5.4, sqrt(area) * 10)),
      text_color = if_else(category %in% c("metabolism", "signaling_kinase", "stress_defense", "hormone", "transport"), "black", "white")
    )

  title_element <- if (is.null(title) || identical(title, "")) {
    element_blank()
  } else {
    element_text(face = "bold", hjust = 0.5, size = 15)
  }

  p <- ggplot(rect_df) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = category), color = "white", linewidth = 1.2) +
    geom_text(aes(x = cx, y = cy, label = label, size = text_size, color = text_color), fontface = "bold", lineheight = 0.9) +
    scale_fill_manual(values = term_colors, breaks = names(term_colors), labels = gsub("_", " ", names(term_colors))) +
    scale_color_identity() +
    scale_size_identity() +
    coord_fixed(expand = FALSE) +
    labs(title = title, fill = NULL) +
    theme_void(base_size = 15) +
    theme(
      plot.title = title_element,
      legend.position = "none"
    )

  ggsave(outfile, p, width = 5.8, height = 5.2, dpi = 300, bg = "white")
}

plot_term_legend <- function(outfile) {
  legend_df <- tibble(
    category = names(term_colors),
    label = gsub("_", " ", names(term_colors)),
    x = rep(1:2, length.out = length(term_colors)),
    y = rep(ceiling(length(term_colors) / 2):1, each = 2)[seq_along(term_colors)]
  )

  p <- ggplot(legend_df, aes(x = x, y = y)) +
    geom_tile(aes(fill = category), width = 0.22, height = 0.22) +
    geom_text(aes(x = x + 0.16, label = label), hjust = 0, size = 5, fontface = "bold") +
    scale_fill_manual(values = term_colors) +
    coord_cartesian(xlim = c(0.85, 2.85), clip = "off") +
    theme_void() +
    theme(legend.position = "none")

  ggsave(outfile, p, width = 7.2, height = 3.8, dpi = 300, bg = "white")
}

keygene_detail <- read_csv(paths$keygene_detail, show_col_types = FALSE)
degkg_neighbors <- read_csv(paths$degKG_neighbors, show_col_types = FALSE)
nondegkg_neighbors <- read_csv(paths$nondegKG_neighbors, show_col_types = FALSE)

keygene_status <- keygene_detail %>%
  count(status = DEG_status, name = "n") %>%
  mutate(status = recode(status, "non-DEG" = "non-DEG", "DEG" = "DEG")) %>%
  arrange(match(status, c("DEG", "non-DEG")))

degkg_status <- degkg_neighbors %>%
  count(status = neighbor_DEG_status, name = "n") %>%
  arrange(match(status, c("DEG", "non-DEG")))

nondegkg_status <- nondegkg_neighbors %>%
  count(status = neighbor_DEG_status, name = "n") %>%
  arrange(match(status, c("DEG", "non-DEG")))

degkg_module_status <- degkg_neighbors %>%
  count(module = neighbor_module, status = neighbor_DEG_status, name = "n") %>%
  complete(module = c("HEAT", "PS", "ROS"), status = c("non-DEG", "DEG"), fill = list(n = 0))

nondegkg_module_status <- nondegkg_neighbors %>%
  count(module = neighbor_module, status = neighbor_DEG_status, name = "n") %>%
  complete(module = c("HEAT", "PS", "ROS"), status = c("non-DEG", "DEG"), fill = list(n = 0))

d_degkg_connected_DEG <- degkg_neighbors %>%
  filter(neighbor_DEG) %>%
  count_keywords()
e_degkg_connected_nonDEG <- degkg_neighbors %>%
  filter(!neighbor_DEG) %>%
  count_keywords()
h_nondegkg_connected_DEG <- nondegkg_neighbors %>%
  filter(neighbor_DEG) %>%
  count_keywords()
i_nondegkg_connected_nonDEG <- nondegkg_neighbors %>%
  filter(!neighbor_DEG) %>%
  count_keywords()

plot_status_pie(
  keygene_status,
  file.path(out_dir, "panel_a_keygene_DEG_status_pie.png"),
  "KeyGene DEG status",
  status_colors_main,
  show_percent = TRUE
)
plot_status_pie(
  degkg_status,
  file.path(out_dir, "panel_b_DEG_keygene_connected_gene_DEG_status_pie.png"),
  "Connected genes: DEG status",
  status_colors_red
)
plot_module_status_bar(
  degkg_module_status,
  file.path(out_dir, "panel_c_DEG_keygene_connected_gene_module_by_DEG_status_bar.png"),
  "Connected genes: module",
  status_colors_red
)
plot_keyword_treemap(
  d_degkg_connected_DEG,
  file.path(out_dir, "panel_d_DEG_keygene_connected_DEG_gene_function_treemap.png"),
  "Connected DEG genes: function"
)
plot_keyword_treemap(
  e_degkg_connected_nonDEG,
  file.path(out_dir, "panel_e_DEG_keygene_connected_nonDEG_gene_function_treemap.png"),
  "Connected non-DEG genes: function"
)
plot_status_pie(
  nondegkg_status,
  file.path(out_dir, "panel_f_nonDEG_keygene_connected_gene_DEG_status_pie.png"),
  "Connected genes: DEG status",
  status_colors_blue
)
plot_module_status_bar(
  nondegkg_module_status,
  file.path(out_dir, "panel_g_nonDEG_keygene_connected_gene_module_by_DEG_status_bar.png"),
  "Connected genes: module",
  status_colors_blue
)
plot_keyword_treemap(
  h_nondegkg_connected_DEG,
  file.path(out_dir, "panel_h_nonDEG_keygene_connected_DEG_gene_function_treemap.png"),
  "Connected DEG genes: function"
)
plot_keyword_treemap(
  i_nondegkg_connected_nonDEG,
  file.path(out_dir, "panel_i_nonDEG_keygene_connected_nonDEG_gene_function_treemap.png"),
  "Connected non-DEG genes: function"
)
plot_term_legend(file.path(out_dir, "panel_function_term_color_legend.png"))

summary_tables <- list(
  panel_a_keygene_status = keygene_status,
  panel_b_DEG_keygene_connected_status = degkg_status,
  panel_c_DEG_keygene_module_status = degkg_module_status,
  panel_d_DEG_keygene_connected_DEG_function = d_degkg_connected_DEG,
  panel_e_DEG_keygene_connected_nonDEG_function = e_degkg_connected_nonDEG,
  panel_f_nonDEG_keygene_connected_status = nondegkg_status,
  panel_g_nonDEG_keygene_module_status = nondegkg_module_status,
  panel_h_nonDEG_keygene_connected_DEG_function = h_nondegkg_connected_DEG,
  panel_i_nonDEG_keygene_connected_nonDEG_function = i_nondegkg_connected_nonDEG
)

for (nm in names(summary_tables)) {
  write_csv(summary_tables[[nm]], file.path(out_dir, paste0(nm, ".csv")))
}

cat("\n===== DONE =====\n")
cat("Output directory:", out_dir, "\n\n")
cat("[a]\n")
print(keygene_status)
cat("\n[b]\n")
print(degkg_status)
cat("\n[f]\n")
print(nondegkg_status)
