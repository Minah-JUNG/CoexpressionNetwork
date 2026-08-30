###############################################################
# DEG-only subcluster module selection from BP ORA terms
###############################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(openxlsx)
})

# === PATHS ===
dir.base <- "/path/to/project"
dir.in   <- file.path(dir.base, "Network_DEGonly/30_Function_Analysis/subcluster")
# ===========================================

file.bp <- file.path(dir.in, "ORA_subcluster_BP_summary.csv")
file.tair <- file.path(dir.in, "subcluster_tair_info.csv")

stopifnot(file.exists(file.bp), file.exists(file.tair))

bp <- read_csv(file.bp, show_col_types = FALSE)
tair_info <- read_csv(file.tair, show_col_types = FALSE)

module_rules <- tibble::tribble(
  ~Module, ~Priority, ~Description, ~Rule,
  "CELL CYCLE", 10,
  "Genes controlling DNA replication, repair, chromosome behavior, and mitotic cell-cycle progression.",
  "dna replication|dna-templated dna replication|cell cycle|mitotic|chromosome|chromatid|cytokinesis|dna damage|dna repair|recombinational repair|g1/s",
  "CELL WALL", 20,
  "Genes directing cell-wall remodeling and secondary wall / phenylpropanoid / lignin-related metabolism.",
  "cell wall|lignin|xylan|hemicellulose|phenylpropanoid|suberin|cutin|wax|glucan|polysaccharide|pectin|secondary metabol",
  "HEAT", 30,
  "Genes mediating heat response, protein folding, chaperone activity, ER stress, and proteostasis.",
  "response to heat|heat acclimation|protein folding|protein refolding|protein maturation|chaperone|endoplasmic reticulum stress|erad pathway|topologically incorrect protein|proteostasis",
  "HORMONE", 40,
  "Genes involved in hormone metabolism and signaling.",
  "hormone|jasmonic|salicylic|ethylene|auxin|gibberellin|abscisic|cytokinin|brassinosteroid|strigolactone",
  "PS", 50,
  "Genes involved in photosynthesis, chloroplast/plastid organization, thylakoid/light reactions, and energy generation.",
  "photosynthesis|photosynthetic|chloroplast|plastid|thylakoid|chlorophyll|light harvesting|electron transport|precursor metabolites and energy|cytochrome complex|pigment metabolic",
  "RNA", 60,
  "Genes governing RNA processing, splicing, ribosome biogenesis, translation, and RNA metabolic control.",
  "rna|mrna|rrna|ribosome|ribonucleoprotein|spliceosome|splicing|translation|translational|nucleolar|rna polymerase|transcription",
  "ROS", 70,
  "Genes linking oxygen/redox stress, hypoxia, detoxification, immune/defense response, and stress signaling.",
  "reactive oxygen|hydrogen peroxide|oxidative stress|hypoxia|decreased oxygen|oxygen levels|detoxification|toxic substance|hypersensitive|immune response|defense response|water deprivation|uv protection",
  "UPS", 80,
  "Genes executing ubiquitin-proteasome and protein catabolic processes.",
  "ubiquitin|proteasom|protein catabolic|polyubiquitination|modification-dependent protein catabolic",
  "VESICLE", 90,
  "Genes mediating vesicular traffic, secretion, Golgi/ER transport, endosome trafficking, and exocytosis.",
  "vesicle|golgi|exocytosis|secretion|secretory|endosome|copii|membrane trafficking|transport, golgi|er to golgi"
)

module_order <- module_rules$Module

assign_modules <- function(description) {
  desc <- str_to_lower(description)
  hits <- module_rules %>%
    rowwise() %>%
    mutate(match = str_detect(desc, Rule)) %>%
    ungroup() %>%
    filter(match) %>%
    arrange(Priority)

  if (nrow(hits) == 0) {
    return(tibble(Module = "OTHER", module_priority = 999, matched_rule = NA_character_))
  }

  tibble(
    Module = hits$Module,
    module_priority = hits$Priority,
    matched_rule = hits$Rule
  )
}

bp_tagged <- bp %>%
  mutate(row_id = row_number()) %>%
  rowwise() %>%
  mutate(module_hit = list(assign_modules(Description))) %>%
  ungroup() %>%
  unnest(module_hit) %>%
  group_by(row_id) %>%
  arrange(module_priority, p.adjust, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    Module = factor(Module, levels = c(module_order, "OTHER")),
    Module = as.character(Module)
  ) %>%
  select(-row_id)

selected_subclusters <- bp_tagged %>%
  filter(Module != "OTHER") %>%
  group_by(subcluster_id, major_cluster) %>%
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
  arrange(factor(Module, levels = module_order), major_cluster, representative_p.adjust)

mapping_sheet <- tair_info %>%
  mutate(
    `function check` = if_else(subcluster_id %in% selected_subclusters$subcluster_id, 1L, 0L)
  ) %>%
  rename(
    `subcluster id` = subcluster_id,
    `total genes` = total_genes,
    `tair genes` = tair_genes
  )

module_summary <- selected_subclusters %>%
  group_by(Module) %>%
  summarise(
    Description = module_rules$Description[match(first(Module), module_rules$Module)],
    `GO terms` = paste(unique(representative_Description), collapse = " | "),
    Pattern = paste(unique(major_cluster), collapse = "; "),
    `Number of subclusters` = n_distinct(subcluster_id),
    `Best p.adjust` = min(representative_p.adjust, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(factor(Module, levels = module_order))

bp_for_excel <- bp_tagged %>%
  select(
    subcluster_id, major_cluster, ID, Description, p.adjust,
    Count, GeneRatio, Module, module_priority, matched_rule
  ) %>%
  arrange(factor(Module, levels = c(module_order, "OTHER")), major_cluster, p.adjust)

output_xlsx <- file.path(
  dir.in,
  "ORA_subcluster_BP_summary_selection_DEGonly_HiSat_2026-05-07.xlsx"
)

write_csv(bp_for_excel, file.path(dir.in, "DEGonly_subcluster_function_module_assignment.csv"))
write_csv(selected_subclusters, file.path(dir.in, "DEGonly_selected_module_subclusters.csv"))
write_csv(module_summary, file.path(dir.in, "DEGonly_module_summary.csv"))
write_csv(module_rules, file.path(dir.in, "DEGonly_module_selection_rules.csv"))

wb <- createWorkbook()
addWorksheet(wb, "subcluster tair mapping")
addWorksheet(wb, "subcluster function")
addWorksheet(wb, "selected module subclusters")
addWorksheet(wb, "module summary")
addWorksheet(wb, "module rules")

writeData(wb, "subcluster tair mapping", mapping_sheet)
writeData(wb, "subcluster function", bp_for_excel)
writeData(wb, "selected module subclusters", selected_subclusters)
writeData(wb, "module summary", module_summary)
writeData(wb, "module rules", module_rules)

header_style <- createStyle(textDecoration = "bold", fgFill = "#D9EAF7", border = "Bottom")
for (sheet in names(wb)) {
  addStyle(wb, sheet, header_style, rows = 1, cols = 1:50, gridExpand = TRUE, stack = TRUE)
  freezePane(wb, sheet, firstRow = TRUE)
  setColWidths(wb, sheet, cols = 1:50, widths = "auto")
}

saveWorkbook(wb, output_xlsx, overwrite = TRUE)

cat("Input BP rows:", nrow(bp), "\n")
cat("BP subclusters:", dplyr::n_distinct(bp$subcluster_id), "\n")
cat("Selected module subclusters:", nrow(selected_subclusters), "\n")
cat("Output xlsx:", output_xlsx, "\n")
cat("\nModule counts:\n")
print(selected_subclusters %>% count(Module, sort = TRUE))
