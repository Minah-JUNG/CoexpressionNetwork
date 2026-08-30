###############################################################
# 03. Module sub-network construction (PS / HEAT / ROS)
###############################################################

rm(list = ls())

library(dplyr)
library(readr)
library(igraph)
library(ggraph)
library(tidyr)
library(stringr)
library(tibble)
library(RColorBrewer)

set.seed(777)

#--------------------------------------------------------------
# 1. Directory
#--------------------------------------------------------------

# === PATHS ===
dir.base <- "/path/to/results"
dir.out  <- file.path(dir.base, "004_Network/module_subnetwork")
# =============================================
if (!dir.exists(dir.out)) dir.create(dir.out, recursive = TRUE)

#--------------------------------------------------------------
# 2. Module definitions
#--------------------------------------------------------------

source("01_module_definition.R")   # axis_tbl, axis_palette, assign_module()

#--------------------------------------------------------------
# 3. Load data
#--------------------------------------------------------------

df_node <- read_csv(
  file.path(dir.base,
            "004_Network/Leiden_Subgrouping_iterative/<DATE>/Iterative_node_membership_thr200.csv")
)

df_cor <- readRDS(
  file.path(dir.base, "003_Correlation/Filtered_correlation_pairs.rds")
)

#--------------------------------------------------------------
# 4. Assign module to each node
#--------------------------------------------------------------

df_axis <- assign_module(df_node) %>% filter(!is.na(axis))

axis_genes <- df_axis %>%
  group_by(axis) %>%
  summarise(n_genes = n_distinct(node))

print(axis_genes)

#--------------------------------------------------------------
# 5. Correlation filtering
#--------------------------------------------------------------

th.cor <- 0.9

df_cor_f <- df_cor %>%
  filter(abs(Correlation) > th.cor)

all_axis_genes <- unique(df_axis$node)

edges_axis <- df_cor_f %>%
  filter(Gene1 %in% all_axis_genes,
         Gene2 %in% all_axis_genes)

#--------------------------------------------------------------
# 6. Build graph
#--------------------------------------------------------------

g <- graph_from_data_frame(edges_axis, directed = FALSE)

axis_map <- df_axis %>%
  distinct(node, axis)

V(g)$axis <- axis_map$axis[match(V(g)$name, axis_map$node)]

#--------------------------------------------------------------
# 7. Remove small components
#--------------------------------------------------------------

comp <- components(g)

comp$membership %>% table %>% as.matrix

# 2025-12-18: 제일 큰 component 만 그림으로
max(comp$csize)

keep_comp <- which(comp$csize > 4)

g <- induced_subgraph(g, V(g)[comp$membership %in% keep_comp])

cat("Graph size:", vcount(g), "nodes,", ecount(g), "edges\n")

#--------------------------------------------------------------
# 8. Centrality metrics
#--------------------------------------------------------------

node_df <- tibble(
  name        = V(g)$name,
  axis        = V(g)$axis,
  degree      = degree(g),
  betweenness = betweenness(g, normalized = TRUE),
  closeness   = closeness(g)
)

cut_deg   <- quantile(node_df$degree, 0.95)
cut_bet   <- quantile(node_df$betweenness, 0.95)
cut_close <- quantile(node_df$closeness, 0.95)

#--------------------------------------------------------------
# 9. Bridge gene detection
#--------------------------------------------------------------

edge_tbl <- igraph::as_data_frame(g, what = "edges") %>%
  mutate(
    axis1 = node_df$axis[match(from, node_df$name)],
    axis2 = node_df$axis[match(to,   node_df$name)]
  )

bridge_genes <- edge_tbl %>%
  filter(axis1 != axis2) %>%
  select(from, to) %>%
  unlist() %>%
  unique()

#--------------------------------------------------------------
# 10. Define KeyGene
#--------------------------------------------------------------

node_df <- node_df %>%
  mutate(
    TopCentrality = degree >= cut_deg |
      betweenness >= cut_bet |
      closeness >= cut_close,
    Bridge = name %in% bridge_genes,
    KeyGene = TopCentrality & Bridge
  )

write_csv(node_df,
          file.path(dir.out, "node_centrality_KeyGene.csv"))

#--------------------------------------------------------------
# 11. Plot
#--------------------------------------------------------------

pal <- axis_palette   # from 000_module_definition.R

V(g)$KeyGene <- node_df$KeyGene[match(V(g)$name, node_df$name)]

layout_fr <- create_layout(g, layout = "fr")

p <- ggraph(layout_fr) +
  geom_edge_link(alpha = 0.15, color = "grey50") +
  geom_node_point(
    aes(fill = axis, shape = as.factor(KeyGene)),
    size = 4, color = "black", stroke = 0.4
  ) +
  geom_node_text(
    aes(label = ifelse(KeyGene, name, "")),
    repel = TRUE, size = 3, fontface = "bold"
  ) +
  scale_fill_manual(values = pal) +
  scale_shape_manual(values = c("FALSE"=21,"TRUE"=24)) +
  theme_graph() +
  labs(title = "Multiple Cluster Network (PS–Heat–ROS)")

ggsave(
  file.path(dir.out, "Multiple_Cluster_network.png"),
  p, width = 12, height = 10, dpi = 300
)

cat("===== ALL DONE =====\n")
