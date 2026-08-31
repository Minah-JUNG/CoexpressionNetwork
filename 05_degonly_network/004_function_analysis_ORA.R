## Cluster/subcluster functional analysis for the DEG-only network.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(openxlsx)
  library(clusterProfiler)
  library(org.At.tair.db)
})

select <- dplyr::select
filter <- dplyr::filter
slice <- dplyr::slice
slice_head <- dplyr::slice_head
arrange <- dplyr::arrange
rename <- dplyr::rename
mutate <- dplyr::mutate
summarise <- dplyr::summarise
first <- dplyr::first

# === PATHS ===
BASE     <- "/path/to/project"
NET      <- file.path(BASE, "Network_DEGonly")
OUT_ROOT <- file.path(NET, "30_Function_Analysis")
BLAST    <- "/path/to/Blastp_TAIR11/blastp_araport11_filtered_best.tsv"
# ===========================================

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

message2 <- function(...) cat(paste0(...), "\n")

load_blast_best <- function() {
  read_tsv(BLAST, show_col_types = FALSE) %>%
    mutate(
      geneID = str_replace(qseqid, "\\.t[0-9]+$", ""),
      tair = str_replace(sseqid, "\\.[0-9]+$", "")
    ) %>%
    arrange(geneID, evalue, desc(bitscore), desc(pident), desc(qcov), desc(scov)) %>%
    group_by(geneID) %>%
    slice(1) %>%
    ungroup() %>%
    select(geneID, tair, qseqid, sseqid, pident, evalue, bitscore, qcov, scov)
}

run_ora_by_group <- function(map_df, group_col, out_dir, prefix) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  tair_info <- map_df %>%
    filter(!is.na(tair)) %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      total_genes = n(),
      tair_genes = n_distinct(tair),
      .groups = "drop"
    ) %>%
    arrange(desc(total_genes))

  write_csv(tair_info, file.path(out_dir, paste0(prefix, "_tair_info.csv")))

  valid_groups <- tair_info %>%
    filter(tair_genes >= 5) %>%
    pull(.data[[group_col]])

  message2("[", prefix, "] valid groups with >=5 TAIR genes: ", length(valid_groups))

  ora_results <- list()
  summaries <- list(BP = list(), CC = list(), MF = list())

  for (i in seq_along(valid_groups)) {
    gid <- valid_groups[i]
    genes <- map_df %>%
      filter(.data[[group_col]] == gid, !is.na(tair)) %>%
      pull(tair) %>%
      unique()

    message2("  [", i, "/", length(valid_groups), "] ", gid, " | TAIR=", length(genes))

    ora <- tryCatch(
      enrichGO(
        gene = genes,
        OrgDb = org.At.tair.db,
        keyType = "TAIR",
        ont = "ALL",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05
      ),
      error = function(e) {
        message2("    ORA failed: ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(ora) && nrow(ora@result) > 0) {
      res <- as.data.frame(ora@result) %>%
        as_tibble() %>%
        mutate(!!group_col := gid, .before = 1)
      ora_results[[as.character(gid)]] <- res

      for (ont in c("BP", "CC", "MF")) {
        top <- res %>%
          filter(ONTOLOGY == ont) %>%
          arrange(p.adjust, pvalue) %>%
          slice_head(n = 10)
        if (nrow(top) > 0) summaries[[ont]][[as.character(gid)]] <- top
      }
    }
  }

  all_ora <- bind_rows(ora_results)
  saveRDS(ora_results, file.path(out_dir, paste0("ORA_", prefix, "_results.rds")))
  write_csv(all_ora, file.path(out_dir, paste0("ORA_", prefix, "_ALL_results.csv")))

  for (ont in c("BP", "CC", "MF")) {
    out <- bind_rows(summaries[[ont]])
    write_csv(out, file.path(out_dir, paste0("ORA_", prefix, "_", ont, "_summary.csv")))
  }

  invisible(list(tair_info = tair_info, ora = all_ora))
}

module_rules <- tibble::tribble(
  ~Module, ~Priority, ~Description, ~Rule,
  "CELL CYCLE", 10, "DNA replication, repair, chromosome behavior, and mitotic cell-cycle progression.",
  "dna replication|dna-templated dna replication|cell cycle|mitotic|chromosome|chromatid|cytokinesis|dna damage|dna repair|recombinational repair|g1/s",
  "CELL WALL", 20, "Cell-wall remodeling and secondary wall / phenylpropanoid / lignin-related metabolism.",
  "cell wall|lignin|xylan|hemicellulose|phenylpropanoid|suberin|cutin|wax|glucan|polysaccharide|pectin|secondary metabol",
  "HEAT", 30, "Heat response, protein folding, chaperone activity, ER stress, and proteostasis.",
  "response to heat|heat acclimation|protein folding|protein refolding|protein maturation|chaperone|endoplasmic reticulum stress|erad pathway|topologically incorrect protein|proteostasis",
  "HORMONE", 40, "Hormone metabolism and signaling.",
  "hormone|jasmonic|salicylic|ethylene|auxin|gibberellin|abscisic|cytokinin|brassinosteroid|strigolactone",
  "PS", 50, "Photosynthesis, chloroplast/plastid organization, thylakoid/light reactions, and energy generation.",
  "photosynthesis|photosynthetic|chloroplast|plastid|thylakoid|chlorophyll|light harvesting|electron transport|precursor metabolites and energy|cytochrome complex|pigment metabolic",
  "RNA", 60, "RNA processing, splicing, ribosome biogenesis, translation, and RNA metabolic control.",
  "rna|mrna|rrna|ribosome|ribonucleoprotein|spliceosome|splicing|translation|translational|nucleolar|rna polymerase|transcription",
  "ROS", 70, "Oxygen/redox stress, hypoxia, detoxification, immune/defense response, and stress signaling.",
  "reactive oxygen|hydrogen peroxide|oxidative stress|hypoxia|decreased oxygen|oxygen levels|detoxification|toxic substance|hypersensitive|immune response|defense response|water deprivation|uv protection",
  "UPS", 80, "Ubiquitin-proteasome and protein catabolic processes.",
  "ubiquitin|proteasom|protein catabolic|polyubiquitination|modification-dependent protein catabolic",
  "VESICLE", 90, "Vesicular traffic, secretion, Golgi/ER transport, endosome trafficking, and exocytosis.",
  "vesicle|golgi|exocytosis|secretion|secretory|endosome|copii|membrane trafficking|transport, golgi|er to golgi"
)

assign_module <- function(description) {
  desc <- str_to_lower(description)
  hit <- module_rules %>%
    rowwise() %>%
    mutate(is_hit = str_detect(desc, Rule)) %>%
    ungroup() %>%
    filter(is_hit) %>%
    arrange(Priority) %>%
    slice_head(n = 1)

  if (nrow(hit) == 0) {
    tibble(Module = "OTHER", module_priority = 999, matched_rule = NA_character_)
  } else {
    tibble(Module = hit$Module, module_priority = hit$Priority, matched_rule = hit$Rule)
  }
}

select_subcluster_modules <- function(sub_out) {
  bp_file <- file.path(sub_out, "ORA_subcluster_BP_summary.csv")
  tair_file <- file.path(sub_out, "subcluster_tair_info.csv")
  if (!file.exists(bp_file) || !file.exists(tair_file)) {
    message2("Module selection skipped: missing BP summary or tair info.")
    return(invisible(NULL))
  }

  bp <- read_csv(bp_file, show_col_types = FALSE)
  tair_info <- read_csv(tair_file, show_col_types = FALSE)
  if (nrow(bp) == 0) {
    message2("Module selection skipped: no BP terms.")
    return(invisible(NULL))
  }

  bp_tagged <- bp %>%
    mutate(row_id = row_number()) %>%
    rowwise() %>%
    mutate(module_hit = list(assign_module(Description))) %>%
    ungroup() %>%
    unnest(module_hit) %>%
    group_by(row_id) %>%
    arrange(module_priority, p.adjust, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    select(-row_id)

  selected <- bp_tagged %>%
    filter(Module != "OTHER") %>%
    group_by(subcluster_id) %>%
    arrange(p.adjust, module_priority, .by_group = TRUE) %>%
    summarise(
      Module = first(Module),
      representative_GO = first(ID),
      representative_Description = first(Description),
      representative_p.adjust = first(p.adjust),
      representative_Count = first(Count),
      representative_GeneRatio = first(GeneRatio),
      matched_terms = paste(unique(Description), collapse = " | "),
      matched_modules = paste(unique(Module), collapse = "; "),
      n_matched_BP_terms = n(),
      .groups = "drop"
    ) %>%
    mutate(major_cluster = str_extract(subcluster_id, "^C[0-9]+")) %>%
    arrange(factor(Module, levels = module_rules$Module), representative_p.adjust)

  module_summary <- selected %>%
    group_by(Module) %>%
    summarise(
      Description = module_rules$Description[match(first(Module), module_rules$Module)],
      `GO terms` = paste(unique(representative_Description), collapse = " | "),
      Pattern = paste(unique(major_cluster), collapse = "; "),
      `Number of subclusters` = n_distinct(subcluster_id),
      `Best p.adjust` = min(representative_p.adjust, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(factor(Module, levels = module_rules$Module))

  mapping <- tair_info %>%
    mutate(`function check` = if_else(subcluster_id %in% selected$subcluster_id, 1L, 0L))

  write_csv(bp_tagged, file.path(sub_out, "DEGonly_subcluster_function_module_assignment.csv"))
  write_csv(selected, file.path(sub_out, "DEGonly_selected_module_subclusters.csv"))
  write_csv(module_summary, file.path(sub_out, "DEGonly_module_summary.csv"))
  write_csv(module_rules, file.path(sub_out, "DEGonly_module_selection_rules.csv"))

  wb <- createWorkbook()
  addWorksheet(wb, "subcluster tair mapping")
  addWorksheet(wb, "subcluster function")
  addWorksheet(wb, "module summary")
  addWorksheet(wb, "selected subclusters")
  addWorksheet(wb, "module rules")
  writeData(wb, "subcluster tair mapping", mapping)
  writeData(wb, "subcluster function", bp_tagged)
  writeData(wb, "module summary", module_summary)
  writeData(wb, "selected subclusters", selected)
  writeData(wb, "module rules", module_rules)
  for (s in names(wb)) {
    freezePane(wb, s, firstRow = TRUE)
    setColWidths(wb, s, cols = 1:50, widths = "auto")
  }
  saveWorkbook(wb, file.path(sub_out, "ORA_subcluster_BP_summary_selection_DEGonly_2026-05-13.xlsx"), overwrite = TRUE)
}

blast_best <- load_blast_best()

cluster_in <- file.path(NET, "10_Cluster/Table_NodeInfo/2026-05-13/Leiden_res1_membership.csv")
cluster_map <- read_csv(cluster_in, show_col_types = FALSE) %>%
  rename(node = Gene, cluster_num = cluster) %>%
  mutate(cluster = paste0("C", cluster_num)) %>%
  left_join(blast_best %>% select(geneID, tair), by = c("node" = "geneID"))

sub_in <- file.path(NET, "20_SubCluster/Leiden_Subgrouping_iterative/2026-05-13/Iterative_node_membership_thr200.csv")
sub_map <- read_csv(sub_in, show_col_types = FALSE) %>%
  mutate(major_cluster = paste0("C", res1_cluster)) %>%
  left_join(blast_best %>% select(geneID, tair), by = c("node" = "geneID"))

message2("========== CLUSTER ORA ==========")
run_ora_by_group(cluster_map, "cluster", file.path(OUT_ROOT, "cluster"), "cluster")

message2("========== SUBCLUSTER ORA ==========")
run_ora_by_group(sub_map, "subcluster_id", file.path(OUT_ROOT, "subcluster"), "subcluster")

# NOTE: Module selection moved to 005_subcluster_module_selection.R.
#       Run that script after this one to assign each subcluster to a
#       module (CELL CYCLE / CELL WALL / HEAT / HORMONE / PS / RNA / ROS /
#       UPS / VESICLE).

message2("DONE: ", OUT_ROOT)
