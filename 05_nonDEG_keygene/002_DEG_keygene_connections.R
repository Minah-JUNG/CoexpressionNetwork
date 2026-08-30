###############################################################
# SET1 DEG KeyGene neighbor analysis
###############################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
  library(tibble)
  library(igraph)
})

# === PATHS ===
out_root <- "/path/to/output_root"
out_dir  <- file.path(out_root,
  paste0("SET1_DEG_KeyGene_connected_nonDEG_analysis_", format(Sys.Date(), "%Y-%m-%d")))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  graph_rds        = "/path/to/Network/Backbone_From_Subclusters/Global_network_backbone_from_subclusters_thr200.rds",
  full_node_table  = "/path/to/Full_vs_DEGonly_subnetwork/full_subnetwork_node_centrality_KeyGene.csv",
  set1_binary      = "/path/to/DEG/DEG_matrix_binary.txt",
  blastp_best      = "/path/to/Blastp_TAIR11/blastp_araport11_filtered_best.tsv",
  araport_function = "/path/to/Araport11_functional_descriptions.txt"
)
# =============================================

comparisons_set1 <- c("T25_H1_T42_H1", "T25_D1_T42_D1")

read_deg_binary <- function(path) {
  comparison_cols <- strsplit(readLines(path, n = 1, warn = FALSE), "\t", fixed = TRUE)[[1]]
  read_tsv(path, skip = 1, col_names = c("gene", comparison_cols), show_col_types = FALSE) %>%
    mutate(across(-gene, ~ suppressWarnings(as.numeric(.x))))
}

clean_text <- function(x) {
  x <- ifelse(is.na(x) | x == "NULL", NA_character_, x)
  iconv(x, from = "", to = "UTF-8", sub = "")
}

sanitize_df <- function(df) {
  df %>%
    mutate(across(where(is.character), ~ iconv(.x, from = "", to = "UTF-8", sub = "")))
}

make_annotation_table <- function() {
  blast <- read_tsv(paths$blastp_best, show_col_types = FALSE) %>%
    mutate(gene = sub("\\.t[0-9]+$", "", qseqid)) %>%
    arrange(gene, evalue, desc(bitscore), desc(qcov), desc(scov)) %>%
    group_by(gene) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      gene,
      araport_id = sseqid,
      blast_pident = pident,
      blast_evalue = evalue,
      blast_bitscore = bitscore,
      blast_qcov = qcov,
      blast_scov = scov
    )

  function_desc <- read_tsv(paths$araport_function, show_col_types = FALSE) %>%
    transmute(
      araport_id = name,
      gene_model_type,
      short_description = clean_text(short_description),
      Curator_summary = clean_text(Curator_summary),
      Computational_description = clean_text(Computational_description)
    )

  blast %>% left_join(function_desc, by = "araport_id")
}

make_keyword_hits <- function(df) {
  patterns <- list(
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

  text <- paste(
    ifelse(is.na(df$short_description), "", df$short_description),
    ifelse(is.na(df$Curator_summary), "", df$Curator_summary),
    ifelse(is.na(df$Computational_description), "", df$Computational_description)
  )

  hit_df <- bind_cols(
    df %>% select(gene),
    as_tibble(lapply(patterns, function(pattern) grepl(pattern, text, ignore.case = TRUE)))
  )

  hit_df
}

summarise_keywords <- function(df) {
  hits <- make_keyword_hits(df)
  hits %>%
    select(-gene) %>%
    summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "category", values_to = "n_gene") %>%
    arrange(desc(n_gene), category)
}

plot_pie <- function(summary_df, outfile, title) {
  plot_df <- summary_df %>%
    mutate(
      pct = n / sum(n) * 100,
      label = paste0(status, "\n", n, " (", sprintf("%.1f", pct), "%)")
    )

  p <- ggplot(plot_df, aes(x = "", y = n, fill = status)) +
    geom_col(width = 1, color = "white", linewidth = 0.4) +
    coord_polar(theta = "y") +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4) +
    scale_fill_manual(values = c("DEG" = "#FDB462", "non-DEG" = "#B3CDE3")) +
    labs(title = title, fill = NULL) +
    theme_void(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")
  ggsave(outfile, p, width = 5.6, height = 5.2, dpi = 300)
}

plot_bar <- function(summary_df, outfile, title, x_col, y_col = "n") {
  p <- ggplot(summary_df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_col(width = 0.72, fill = "#B3CDE3", color = "white", linewidth = 0.3) +
    labs(title = title, x = NULL, y = "Gene count") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
  ggsave(outfile, p, width = 6.4, height = 4.6, dpi = 300)
}

get_neighbor_pairs <- function(sg, keygenes) {
  keygenes <- intersect(keygenes, V(sg)$name)
  bind_rows(lapply(keygenes, function(keygene) {
    nb <- names(neighbors(sg, keygene, mode = "all"))
    tibble(DEG_keygene = keygene, neighbor_gene = nb)
  }))
}

node_table <- read_csv(paths$full_node_table, show_col_types = FALSE)
annotation <- make_annotation_table()

set1_binary <- read_deg_binary(paths$set1_binary)
set1_genes <- set1_binary %>%
  filter(rowSums(across(all_of(comparisons_set1)), na.rm = TRUE) > 0) %>%
  pull(gene)

cat("Loading graph...\n")
g <- readRDS(paths$graph_rds)
subgraph_nodes <- intersect(node_table$name, V(g)$name)
sg <- induced_subgraph(g, vids = subgraph_nodes)
cat("Subgraph:", vcount(sg), "nodes,", ecount(sg), "edges\n")

deg_keygenes <- node_table %>%
  filter(KeyGene, name %in% set1_genes) %>%
  transmute(
    gene = name,
    module = axis,
    degree,
    betweenness,
    closeness,
    TopCentrality,
    Bridge,
    DEG = TRUE
  ) %>%
  left_join(annotation, by = "gene")

neighbor_pairs <- get_neighbor_pairs(sg, deg_keygenes$gene)

vertex_info <- as_data_frame(sg, what = "vertices") %>%
  as_tibble() %>%
  transmute(
    gene = name,
    graph_cluster = cluster,
    graph_Degree = Degree,
    graph_KeyGene = KeyGene,
    original_cluster,
    subcluster_id,
    recovered_cluster_id,
    recovered_cluster
  )

neighbor_edges <- neighbor_pairs %>%
  left_join(deg_keygenes %>% select(DEG_keygene = gene, keygene_module = module), by = "DEG_keygene") %>%
  left_join(node_table %>% select(neighbor_gene = name, neighbor_module = axis, neighbor_degree = degree, neighbor_KeyGene = KeyGene), by = "neighbor_gene") %>%
  left_join(vertex_info %>% rename(neighbor_gene = gene), by = "neighbor_gene") %>%
  mutate(
    neighbor_DEG = neighbor_gene %in% set1_genes,
    neighbor_DEG_status = if_else(neighbor_DEG, "DEG", "non-DEG")
  ) %>%
  left_join(annotation %>% rename(neighbor_gene = gene), by = "neighbor_gene")

unique_neighbors <- neighbor_edges %>%
  arrange(neighbor_gene, desc(neighbor_DEG), neighbor_module) %>%
  group_by(neighbor_gene) %>%
  summarise(
    neighbor_module = first(neighbor_module),
    neighbor_DEG = first(neighbor_DEG),
    neighbor_DEG_status = first(neighbor_DEG_status),
    n_DEG_keygene_connected = n_distinct(DEG_keygene),
    connected_DEG_keygenes = paste(sort(unique(DEG_keygene)), collapse = ";"),
    graph_cluster = first(graph_cluster),
    original_cluster = first(original_cluster),
    subcluster_id = first(subcluster_id),
    recovered_cluster_id = first(recovered_cluster_id),
    recovered_cluster = first(recovered_cluster),
    araport_id = first(araport_id),
    blast_pident = first(blast_pident),
    blast_evalue = first(blast_evalue),
    blast_bitscore = first(blast_bitscore),
    blast_qcov = first(blast_qcov),
    blast_scov = first(blast_scov),
    gene_model_type = first(gene_model_type),
    short_description = first(short_description),
    Curator_summary = first(Curator_summary),
    Computational_description = first(Computational_description),
    .groups = "drop"
  )

connected_nonDEG <- unique_neighbors %>%
  filter(!neighbor_DEG) %>%
  rename(gene = neighbor_gene, module = neighbor_module)

neighbor_status_summary <- unique_neighbors %>%
  count(status = neighbor_DEG_status, name = "n") %>%
  arrange(status)

connection_status_summary <- neighbor_edges %>%
  count(status = neighbor_DEG_status, name = "n_connection") %>%
  arrange(status)

connected_nonDEG_module_summary <- connected_nonDEG %>%
  count(module, name = "n") %>%
  arrange(desc(n))

connected_nonDEG_keyword_summary <- summarise_keywords(connected_nonDEG)

connected_nonDEG_keyword_by_module <- make_keyword_hits(connected_nonDEG) %>%
  left_join(connected_nonDEG %>% select(gene, module), by = "gene") %>%
  pivot_longer(
    cols = -c(gene, module),
    names_to = "category",
    values_to = "hit"
  ) %>%
  filter(hit) %>%
  count(module, category, name = "n_gene") %>%
  arrange(module, desc(n_gene), category)

deg_keygene_neighbor_summary <- neighbor_edges %>%
  group_by(DEG_keygene, keygene_module) %>%
  summarise(
    n_neighbor = n_distinct(neighbor_gene),
    n_DEG_neighbor = n_distinct(neighbor_gene[neighbor_DEG]),
    n_nonDEG_neighbor = n_distinct(neighbor_gene[!neighbor_DEG]),
    pct_DEG_neighbor = n_DEG_neighbor / n_neighbor * 100,
    n_HEAT_neighbor = n_distinct(neighbor_gene[neighbor_module == "HEAT"]),
    n_PS_neighbor = n_distinct(neighbor_gene[neighbor_module == "PS"]),
    n_ROS_neighbor = n_distinct(neighbor_gene[neighbor_module == "ROS"]),
    .groups = "drop"
  ) %>%
  left_join(deg_keygenes %>% select(DEG_keygene = gene, araport_id, short_description, Curator_summary), by = "DEG_keygene") %>%
  arrange(desc(n_nonDEG_neighbor), desc(n_neighbor), DEG_keygene)

overall_summary <- tibble(
  set_id = "SET1",
  n_DEG_keygene = nrow(deg_keygenes),
  n_neighbor_edge = nrow(neighbor_edges),
  n_unique_neighbor = nrow(unique_neighbors),
  n_unique_neighbor_DEG = sum(unique_neighbors$neighbor_DEG),
  n_unique_neighbor_nonDEG = sum(!unique_neighbors$neighbor_DEG),
  pct_unique_neighbor_DEG = n_unique_neighbor_DEG / n_unique_neighbor * 100,
  n_connection_to_DEG = sum(neighbor_edges$neighbor_DEG),
  n_connection_to_nonDEG = sum(!neighbor_edges$neighbor_DEG),
  pct_connection_to_DEG = n_connection_to_DEG / n_neighbor_edge * 100
)

plot_pie(
  neighbor_status_summary,
  file.path(out_dir, "SET1_DEG_keygene_unique_neighbor_DEG_status_pie.png"),
  "Genes connected to DEG key genes: DEG status"
)

plot_bar(
  connected_nonDEG_module_summary,
  file.path(out_dir, "SET1_DEG_keygene_connected_nonDEG_module_bar.png"),
  "non-DEG genes connected to DEG key genes: module",
  "module"
)

write_csv(sanitize_df(overall_summary), file.path(out_dir, "SET1_DEG_keygene_neighbor_overall_summary.csv"))
write_csv(sanitize_df(deg_keygenes), file.path(out_dir, "SET1_DEG_keygene_function_annotation.csv"))
write_csv(sanitize_df(neighbor_edges), file.path(out_dir, "SET1_DEG_keygene_neighbor_edges_detail.csv"))
write_csv(sanitize_df(unique_neighbors), file.path(out_dir, "SET1_DEG_keygene_unique_connected_genes_annotation.csv"))
write_csv(sanitize_df(connected_nonDEG), file.path(out_dir, "SET1_DEG_keygene_unique_connected_nonDEG_genes_function_annotation.csv"))
write_csv(sanitize_df(neighbor_status_summary), file.path(out_dir, "SET1_DEG_keygene_unique_neighbor_DEG_status_summary.csv"))
write_csv(sanitize_df(connection_status_summary), file.path(out_dir, "SET1_DEG_keygene_connection_DEG_status_summary.csv"))
write_csv(sanitize_df(connected_nonDEG_module_summary), file.path(out_dir, "SET1_DEG_keygene_connected_nonDEG_module_summary.csv"))
write_csv(sanitize_df(connected_nonDEG_keyword_summary), file.path(out_dir, "SET1_DEG_keygene_connected_nonDEG_function_keyword_summary.csv"))
write_csv(sanitize_df(connected_nonDEG_keyword_by_module), file.path(out_dir, "SET1_DEG_keygene_connected_nonDEG_function_keyword_by_module.csv"))
write_csv(sanitize_df(deg_keygene_neighbor_summary), file.path(out_dir, "SET1_DEG_keygene_neighbor_summary.csv"))

wb <- createWorkbook()
sheets <- list(
  overall_summary = overall_summary,
  DEG_keygenes = deg_keygenes,
  unique_neighbors = unique_neighbors,
  connected_nonDEG = connected_nonDEG,
  unique_status_summary = neighbor_status_summary,
  connection_status_summary = connection_status_summary,
  nonDEG_module_summary = connected_nonDEG_module_summary,
  nonDEG_function_keywords = connected_nonDEG_keyword_summary,
  nonDEG_keywords_by_module = connected_nonDEG_keyword_by_module,
  DEG_keygene_summary = deg_keygene_neighbor_summary
)

for (nm in names(sheets)) {
  addWorksheet(wb, substr(nm, 1, 31))
  writeData(wb, substr(nm, 1, 31), sanitize_df(sheets[[nm]]))
}

for (sheet in names(wb)) {
  freezePane(wb, sheet, firstRow = TRUE)
  setColWidths(wb, sheet, cols = 1:80, widths = "auto")
}

saveWorkbook(wb, file.path(out_dir, "SET1_DEG_KeyGene_connected_nonDEG_analysis.xlsx"), overwrite = TRUE)

cat("\n===== DONE =====\n")
cat("Output directory:", out_dir, "\n\n")
print(overall_summary)
cat("\n[Unique connected genes: DEG status]\n")
print(neighbor_status_summary)
cat("\n[Connected non-DEG genes: module]\n")
print(connected_nonDEG_module_summary)
cat("\n[Connected non-DEG genes: functional keyword]\n")
print(connected_nonDEG_keyword_summary)
