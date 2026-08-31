###############################################################
# DEG-only network visualization
# 12 samples | Leiden res=1 | |cor| > 0.95
###############################################################

rm(list = ls())
set.seed(12345) # node 위치 고정

library(dplyr)
library(readr)
library(igraph)
library(ggraph)
library(tidyr)
library(stringr)
library(tibble)
library(tidygraph)
library(viridis)

samples <- 12

# === PATHS ===
NET      <- "/path/to/Network_DEGonly"
OUT_ROOT <- "/path/to/output_root"
# ===========================================

#--------------------------------------------------------------
# Directories
#--------------------------------------------------------------
dir.out <- file.path(OUT_ROOT, "10_Cluster/Network_Figure")
if (!dir.exists(dir.out)) dir.create(dir.out, recursive = TRUE)

#--------------------------------------------------------------
# Load files
#--------------------------------------------------------------
file_node <- file.path(NET, "10_Cluster/Table_NodeInfo/<DATE>/Node_annotation_Leiden_res1.csv")
file_corr <- file.path(NET, "Filtered_correlation_pairs.rds")

node_df <- read_csv(file_node)   # name, cluster, degree, betweenness, closeness, HubGene, ...
corr_df <- readRDS(file_corr)    # Gene_1, Gene_2, Correlation, Pvalue (스크립트1 setNames 기준)

#--------------------------------------------------------------
# Step 1. Prepare node metadata
#--------------------------------------------------------------
major_clusters <- node_df %>%
  count(cluster, name = "node_count") %>%
  filter(node_count >= 10) %>%
  pull(cluster)

node_df <- node_df %>%
  mutate(
    major_cluster = if_else(cluster %in% major_clusters,
                            paste0("C_", cluster),
                            "Other")
  ) %>%
  rename(name = name)   # 이미 'name' 컬럼이므로 유지

#--------------------------------------------------------------
# Step 2. Filter correlation by threshold
#--------------------------------------------------------------
th.cor <- 0.95

corr_ft <- corr_df %>%
  filter(abs(Correlation) >= th.cor)

cat("Pairs after |cor| >=", th.cor, ":", nrow(corr_ft), "\n")

# degree=1 노드(양쪽 모두 singleton) 제거
deg_table   <- table(c(corr_ft$Gene1, corr_ft$Gene2))
single_nodes <- names(deg_table[deg_table == 1])

corr_ft <- corr_ft %>%
  filter(!(Gene1 %in% single_nodes & Gene2 %in% single_nodes))

cat("Pairs after singleton removal:", nrow(corr_ft), "\n")

#--------------------------------------------------------------
# Step 3. Build graph & join node metadata
#--------------------------------------------------------------
g <- graph_from_data_frame(corr_ft[, c("Gene1", "Gene2")], directed = FALSE) %>%
  simplify(remove.multiple = TRUE, remove.loops = TRUE)

tg <- g %>%
  as_tbl_graph() %>%
  activate(nodes) %>%
  left_join(node_df, by = "name")

#--------------------------------------------------------------
# Step 4. Remove small connected components (keep component size > 10)
#--------------------------------------------------------------
tg <- tg %>%
  activate(nodes) %>%
  mutate(component = group_components()) %>%
  group_by(component) %>%
  filter(n() > 10) %>%
  ungroup()

node_count <- tg %>% activate(nodes) %>% as_tibble() %>% nrow()
edge_count <- tg %>% activate(edges) %>% as_tibble() %>% nrow()
cat("Final graph — nodes:", node_count, "| edges:", edge_count, "\n")

#--------------------------------------------------------------
# Step 5. Color palette
#--------------------------------------------------------------
cluster_levels <- tg %>%
  activate(nodes) %>%
  pull(major_cluster) %>%
  unique() %>%
  sort()

major_only <- cluster_levels[cluster_levels != "Other"]
n_major    <- length(major_only)

pal_colors <- setNames(
  c(viridis(n_major, option = "viridis"), "grey75"),
  c(major_only, "Other")
)

#--------------------------------------------------------------
# Step 6. Layout (Fruchterman-Reingold)
#--------------------------------------------------------------
set.seed(12345)
lay <- create_layout(tg, layout = "fr", niter = 500)

lay$x <- lay$x * 1.5
lay$y <- lay$y * 1.5

#--------------------------------------------------------------
# Step 7. Draw
#--------------------------------------------------------------
p <- ggraph(lay) +
  geom_edge_link(alpha = 0.12, color = "grey60", linewidth = 0.2) +
  geom_node_point(
    aes(fill = major_cluster),
    shape  = 21,
    size   = 1.8,
    color  = "black",
    stroke = 0.3,
    alpha  = 0.75
  ) +
  scale_fill_manual(
    values = pal_colors,
    name   = "Cluster",
    breaks = c(major_only, "Other")   # Other를 범례 마지막으로
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 5))) +
  theme_graph(base_family = "sans") +
  labs(title = paste0("Global Co-expression Network  (|r| ≥ ", th.cor, ", ", samples, " samples)"))

#--------------------------------------------------------------
# Step 8. Save
#--------------------------------------------------------------
ggsave(
  file.path(dir.out, paste0("Network_global_cor", th.cor, ".png")),
  p, width = 12, height = 10, dpi = 300
)

ggsave(
  file.path(dir.out, paste0("Network_global_cor", th.cor, ".pdf")),
  p, width = 12, height = 10
)

cat("Saved to:", dir.out, "\n")
