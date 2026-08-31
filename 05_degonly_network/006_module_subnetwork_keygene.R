###############################################################
# Full-network vs DEG-only HEAT/ROS module-based KeyGene analysis
###############################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(igraph)
  library(stringr)
  library(tibble)
  library(openxlsx)
})

# === PATHS ===
base    <- "/path/to/project"
cmp_dir <- file.path(base, "NetworkComparison_Full_vs_DEGonly")
net_dir <- file.path(base, "Network_DEGonly")
# ===========================================
dir.create(cmp_dir, recursive = TRUE, showWarnings = FALSE)

th_cor_full    <- 0.9
target_modules <- c("HEAT", "ROS")
interest_genes <- c("g13850", "g35527", "g26591", "g42276", "g23822")

paths <- list(
  full_keygene_csv        = "/path/to/full_network/KeyGene_list/node_centrality_KeyGene.csv",
  full_membership_csv     = "/path/to/full_network/Iterative_node_membership_thr200.csv",
  full_cor_rds            = "/path/to/full_network/Filtered_correlation_pairs.rds",
  deg_graph_rds           = file.path(net_dir, "10_Cluster/Leiden_Results/<DATE>/g_tg_Leiden_FINAL_res1.rds"),
  deg_membership_csv      = file.path(net_dir, "20_SubCluster/Leiden_Subgrouping_iterative/<DATE>/Iterative_node_membership_thr200.csv"),
  deg_selected_module_csv = file.path(net_dir, "30_Function_Analysis/subcluster/DEGonly_selected_module_subclusters.csv"),
  blast_tsv               = "/path/to/Blastp_TAIR11/blastp_araport11_filtered_best.tsv",
  araport_desc            = "/path/to/Araport11_functional_descriptions.txt"
)

standardize_edges <- function(edge_df) {
  edge_df %>%
    rename_with(~ "Gene1", matches("^Gene1$|^Gene_1$|^from$", ignore.case = TRUE)) %>%
    rename_with(~ "Gene2", matches("^Gene2$|^Gene_2$|^to$", ignore.case = TRUE)) %>%
    rename_with(~ "Correlation", matches("^Correlation$|^cor$", ignore.case = TRUE)) %>%
    select(Gene1, Gene2, Correlation, everything()) %>%
    mutate(
      Gene1 = as.character(Gene1),
      Gene2 = as.character(Gene2)
    )
}

build_keygene_table <- function(graph, label) {
  comp <- components(graph)
  keep_comp <- which(comp$csize > 4)
  graph <- induced_subgraph(graph, V(graph)[comp$membership %in% keep_comp])

  node_df <- tibble(
    name = V(graph)$name,
    module = V(graph)$module,
    subcluster_id = V(graph)$subcluster_id,
    degree = as.integer(degree(graph)),
    betweenness = betweenness(graph, normalized = TRUE),
    closeness = closeness(graph)
  )

  cut_deg <- quantile(node_df$degree, 0.95, na.rm = TRUE)
  cut_bet <- quantile(node_df$betweenness, 0.95, na.rm = TRUE)
  cut_close <- quantile(node_df$closeness, 0.95, na.rm = TRUE)

  edge_tbl <- igraph::as_data_frame(graph, what = "edges") %>%
    as_tibble() %>%
    mutate(
      module1 = node_df$module[match(from, node_df$name)],
      module2 = node_df$module[match(to, node_df$name)]
    )

  bridge_genes <- edge_tbl %>%
    filter(module1 != module2) %>%
    select(from, to) %>%
    unlist(use.names = FALSE) %>%
    unique()

  node_df <- node_df %>%
    mutate(
      TopCentrality = degree >= cut_deg |
        betweenness >= cut_bet |
        closeness >= cut_close,
      Bridge = name %in% bridge_genes,
      KeyGene = TopCentrality & Bridge
    ) %>%
    arrange(module, desc(KeyGene), desc(degree), name)

  summary_df <- tibble(
    network = label,
    modules = paste(target_modules, collapse = "/"),
    n_nodes = vcount(graph),
    n_edges = ecount(graph),
    n_HEAT = sum(node_df$module == "HEAT", na.rm = TRUE),
    n_ROS = sum(node_df$module == "ROS", na.rm = TRUE),
    n_bridge_genes = sum(node_df$Bridge),
    n_topcentrality = sum(node_df$TopCentrality),
    n_KeyGene = sum(node_df$KeyGene),
    cutoff_degree_95pct = as.numeric(cut_deg),
    cutoff_betweenness_95pct = as.numeric(cut_bet),
    cutoff_closeness_95pct = as.numeric(cut_close)
  )

  list(graph = graph, node_df = node_df, edge_tbl = edge_tbl, summary = summary_df)
}

build_full_heat_ros <- function() {
  axis_tbl <- tribble(
    ~module, ~pattern,
    "HEAT", "^C_6_2_|^C_6_3|^C_6_4|^C_1_2_2_",
    "ROS",  "^C_1_3_1_1_1|^C_1_3_2_"
  )

  membership <- read_csv(paths$full_membership_csv, show_col_types = FALSE) %>%
    mutate(
      module = case_when(
        str_detect(subcluster_id, axis_tbl$pattern[axis_tbl$module == "HEAT"]) ~ "HEAT",
        str_detect(subcluster_id, axis_tbl$pattern[axis_tbl$module == "ROS"]) ~ "ROS",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(module)) %>%
    distinct(name = node, module, subcluster_id)

  cor_df <- readRDS(paths$full_cor_rds) %>%
    as_tibble() %>%
    standardize_edges() %>%
    filter(abs(Correlation) > th_cor_full)

  edges <- cor_df %>%
    semi_join(membership, by = c("Gene1" = "name")) %>%
    semi_join(membership, by = c("Gene2" = "name")) %>%
    distinct(Gene1, Gene2, .keep_all = TRUE)

  graph <- graph_from_data_frame(
    edges %>% select(from = Gene1, to = Gene2, Correlation),
    directed = FALSE,
    vertices = membership %>% select(name, module, subcluster_id)
  )

  build_keygene_table(graph, "Full_network_HEAT_ROS")
}

build_degonly_heat_ros <- function() {
  graph_obj <- readRDS(paths$deg_graph_rds)
  graph <- if (inherits(graph_obj, "igraph")) graph_obj else as.igraph(graph_obj)

  membership <- read_csv(paths$deg_membership_csv, show_col_types = FALSE) %>%
    mutate(subcluster_id = as.character(subcluster_id))

  selected_modules <- read_csv(paths$deg_selected_module_csv, show_col_types = FALSE) %>%
    filter(Module %in% target_modules) %>%
    select(subcluster_id, Module) %>%
    distinct()

  vertices <- membership %>%
    inner_join(selected_modules, by = "subcluster_id") %>%
    transmute(name = as.character(node), module = Module, subcluster_id) %>%
    distinct()

  graph <- induced_subgraph(graph, intersect(V(graph)$name, vertices$name))
  V(graph)$module <- vertices$module[match(V(graph)$name, vertices$name)]
  V(graph)$subcluster_id <- vertices$subcluster_id[match(V(graph)$name, vertices$name)]

  res <- build_keygene_table(graph, "DEGonly_network_HEAT_ROS")
  res$selected_modules <- selected_modules
  res
}

load_annotation <- function() {
  clean_null <- function(x) {
    x <- as.character(x)
    x[x %in% c("", "NULL", "NA", "N/A")] <- NA_character_
    x
  }

  blast <- read_tsv(paths$blast_tsv, show_col_types = FALSE) %>%
    mutate(
      gene = str_replace(qseqid, "\\.t[0-9]+$", ""),
      araport_model = sseqid,
      araport_locus = str_replace(sseqid, "\\.[0-9]+$", "")
    ) %>%
    arrange(gene, evalue, desc(bitscore), desc(pident), desc(qcov), desc(scov)) %>%
    group_by(gene) %>%
    slice(1) %>%
    ungroup() %>%
    select(name = gene, blast_query = qseqid, araport_model, araport_locus, pident, align_len, evalue, bitscore, qcov, scov)

  araport <- read_tsv(paths$araport_desc, show_col_types = FALSE) %>%
    mutate(
      araport_model = name,
      araport_locus = str_replace(name, "\\.[0-9]+$", ""),
      Curator_summary = clean_null(Curator_summary)
    )

  by_model <- araport %>% select(araport_model, Curator_summary_model = Curator_summary)
  by_locus <- araport %>%
    filter(!is.na(Curator_summary)) %>%
    group_by(araport_locus) %>%
    summarise(Curator_summary_locus = first(Curator_summary), .groups = "drop")

  blast %>%
    left_join(by_model, by = "araport_model") %>%
    left_join(by_locus, by = "araport_locus") %>%
    mutate(Curator_summary = coalesce(Curator_summary_model, Curator_summary_locus)) %>%
    select(name, blast_query, araport_model, araport_locus, pident, align_len, evalue, bitscore, qcov, scov, Curator_summary)
}

yes_no <- function(x) ifelse(!is.na(x) & x, "Keygene", "-")
fmt_int <- function(x) format(as.numeric(x), big.mark = ",", scientific = FALSE, trim = TRUE)
sanitize_df <- function(df) {
  df %>%
    mutate(across(where(is.character), ~ iconv(.x, from = "", to = "UTF-8", sub = "")))
}

full_res <- build_full_heat_ros()
deg_res <- build_degonly_heat_ros()
annotation <- load_annotation()

# Full-network module key genes remain the original PS/HEAT/ROS result.
# DEG-only is recalculated using HEAT/ROS only, because PS was not selected

full_nodes <- read_csv(paths$full_keygene_csv, show_col_types = FALSE) %>%
  transmute(
    name,
    full_module = axis,
    full_degree = as.numeric(degree),
    full_betweenness = as.numeric(betweenness),
    full_closeness = as.numeric(closeness),
    full_TopCentrality = as.logical(TopCentrality),
    full_Bridge = as.logical(Bridge),
    full_KeyGene = as.logical(KeyGene)
  )

deg_nodes <- deg_res$node_df %>%
  rename(deg_module = module, deg_degree = degree, deg_betweenness = betweenness,
         deg_closeness = closeness, deg_TopCentrality = TopCentrality,
         deg_Bridge = Bridge, deg_KeyGene = KeyGene)

full_key <- full_nodes %>% filter(full_KeyGene)
deg_key <- deg_nodes %>% filter(deg_KeyGene)

module_venn <- tibble(
  comparison = "Key genes for module-based networks (Full HEAT/ROS/PS; DEG-only HEAT/ROS)",
  full_network_only = length(setdiff(full_key$name, deg_key$name)),
  shared = length(intersect(full_key$name, deg_key$name)),
  deg_only_network_only = length(setdiff(deg_key$name, full_key$name)),
  full_network_total = full_network_only + shared,
  deg_only_network_total = deg_only_network_only + shared
)

module_union <- full_join(
  full_key %>%
    select(name, full_module, full_degree, full_betweenness, full_closeness,
           full_TopCentrality, full_Bridge, full_KeyGene),
  deg_key %>%
    select(name, deg_module, deg_subcluster_id = subcluster_id, deg_degree,
           deg_betweenness, deg_closeness, deg_TopCentrality, deg_Bridge, deg_KeyGene),
  by = "name"
) %>%
  mutate(
    module_overlap_category = case_when(
      !is.na(full_KeyGene) & !is.na(deg_KeyGene) ~ "Shared_module_KeyGene",
      !is.na(full_KeyGene) & is.na(deg_KeyGene) ~ "Full_module_only_KeyGene",
      is.na(full_KeyGene) & !is.na(deg_KeyGene) ~ "DEGonly_module_only_KeyGene",
      TRUE ~ "Not_module_KeyGene"
    ),
    .after = name
  ) %>%
  left_join(annotation, by = "name") %>%
  arrange(
    factor(module_overlap_category, levels = c("Shared_module_KeyGene", "Full_module_only_KeyGene", "DEGonly_module_only_KeyGene")),
    name
  ) %>%
  sanitize_df()

deg_module_degree <- setNames(deg_res$node_df$degree, deg_res$node_df$name)
candidate_table <- tibble(Gene = interest_genes) %>%
  left_join(
    full_nodes %>%
      filter(name %in% interest_genes) %>%
      transmute(Gene = name, full_module, full_KeyGene),
    by = "Gene"
  ) %>%
  left_join(
    deg_nodes %>%
      filter(name %in% interest_genes) %>%
      transmute(
        Gene = name,
        deg_module,
        deg_degree,
        deg_TopCentrality,
        deg_Bridge,
        deg_KeyGene
      ),
    by = "Gene"
  ) %>%
  mutate(
    `Full-network` = yes_no(full_KeyGene),
    `DEGonly-network` = yes_no(deg_KeyGene),
    `DEGonly module` = coalesce(deg_module, "-"),
    `Connected gene` = if_else(is.na(deg_degree), "-", fmt_int(deg_degree)),
    Notes = case_when(
      is.na(deg_degree) ~ "Not in selected HEAT/ROS module",
      deg_KeyGene ~ "",
      !deg_TopCentrality & !deg_Bridge ~ "Low centrality; no HEAT-ROS bridge",
      !deg_TopCentrality ~ "Low centrality",
      !deg_Bridge ~ "No HEAT-ROS bridge",
      TRUE ~ ""
    )
  ) %>%
  select(Gene, `Full-network`, `DEGonly-network`, `DEGonly module`, `Connected gene`, Notes)

candidate_table <- sanitize_df(candidate_table)

summary_out <- sanitize_df(bind_rows(full_res$summary, deg_res$summary))
selected_count_out <- sanitize_df(deg_res$selected_modules %>% count(Module, name = "n_subclusters"))
module_venn <- sanitize_df(module_venn)

write_csv(sanitize_df(full_res$node_df), file.path(cmp_dir, "Node_centrality_KeyGene_HEAT_ROS_Full_network.csv"))
write_csv(sanitize_df(full_nodes), file.path(cmp_dir, "Node_centrality_KeyGene_reference_Full_HEAT_ROS_PS_network.csv"))
write_csv(sanitize_df(deg_res$node_df), file.path(cmp_dir, "Node_centrality_KeyGene_HEAT_ROS_DEGonly_network.csv"))
write_csv(summary_out, file.path(cmp_dir, "Summary_KeyGene_HEAT_ROS_module_network.csv"))
write_csv(selected_count_out, file.path(cmp_dir, "DEGonly_HEAT_ROS_selected_subcluster_count.csv"))
write_csv(module_venn, file.path(cmp_dir, "Module_HEAT_ROS_keygene_venn_counts.csv"))
write_csv(module_union, file.path(cmp_dir, "Module_HEAT_ROS_KeyGene_comparison_with_Curator_summary.csv"))
write_csv(candidate_table, file.path(cmp_dir, "Candidate_keygene_comparison_HEAT_ROS_module.csv"))

wb <- createWorkbook()
tables <- list(
  HEAT_ROS_summary = summary_out,
  Venn_counts = module_venn,
  Module_KeyGene_union = module_union,
  Full_module_only_KeyGenes = filter(module_union, module_overlap_category == "Full_module_only_KeyGene"),
  DEGonly_module_only_KeyGenes = filter(module_union, module_overlap_category == "DEGonly_module_only_KeyGene"),
  Shared_module_KeyGenes = filter(module_union, module_overlap_category == "Shared_module_KeyGene"),
  Candidate_table = candidate_table,
  DEGonly_selected_subclusters = sanitize_df(deg_res$selected_modules)
)

for (nm in names(tables)) {
  addWorksheet(wb, nm)
  writeData(wb, nm, tables[[nm]])
}
saveWorkbook(wb, file.path(cmp_dir, "Module_HEAT_ROS_KeyGene_comparison_with_Curator_summary.xlsx"), overwrite = TRUE)

cat("\n[HEAT/ROS module summary]\n")
print(bind_rows(full_res$summary, deg_res$summary))
cat("\n[HEAT/ROS module keygene Venn]\n")
print(module_venn)
cat("\n[Candidate genes]\n")
print(candidate_table)
cat("\n===== ALL DONE =====\n")
