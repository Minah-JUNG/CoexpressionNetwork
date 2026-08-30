###############################################################
# 004. Inter-module mean correlation heatmap
###############################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# === PATHS ===
dir.out   <- "/path/to/results/004_Network/module_subnetwork"
file.cor  <- "/path/to/results/003_Correlation/Filtered_correlation_pairs.rds"
file.meta <- file.path(dir.out, "node_centrality_KeyGene.csv")
# =============================================
if (!dir.exists(dir.out)) dir.create(dir.out, recursive = TRUE)

# Load
cor_pairs <- readRDS(file.cor)
meta      <- read.csv(file.meta)

# gene -> axis map (named vector)
gene_axis <- meta %>% select(name, axis) %>% deframe()

# Annotate each pair with its endpoint axes, drop pairs outside modules
cor_pairs <- cor_pairs %>%
  mutate(
    axis1 = gene_axis[Gene1],
    axis2 = gene_axis[Gene2]
  ) %>%
  filter(!is.na(axis1), !is.na(axis2))

axes <- c("HEAT", "PS", "ROS")

# Per axis pair: mean |Pearson r|
mat_df <- expand.grid(A = axes, B = axes) %>%
  rowwise() %>%
  mutate(
    mean_cor = {
      if (A == B) {
        vals <- cor_pairs %>%
          filter(axis1 == A & axis2 == A) %>%
          pull(Correlation)
      } else {
        vals <- cor_pairs %>%
          filter((axis1 == A & axis2 == B) | (axis1 == B & axis2 == A)) %>%
          pull(Correlation)
      }
      if (length(vals) == 0) NA else mean(abs(vals), na.rm = TRUE)
    }
  ) %>%
  ungroup()

mat_df$A <- factor(mat_df$A, levels = axes)
mat_df$B <- factor(mat_df$B, levels = rev(axes))

# Heatmap
p <- ggplot(mat_df, aes(x = A, y = B, fill = mean_cor)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(mean_cor, 3)), size = 4.5, fontface = "bold") +
  scale_fill_gradient(
    low = "#FEE5D9", high = "#CC2222",
    limits = c(0.77, 0.83),
    name = "|r|"
  ) +
  labs(x = NULL, y = NULL, title = "Inter-module mean correlation coefficient") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text   = element_text(size = 12, color = "black"),
    panel.grid  = element_blank(),
    plot.title  = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.position = "right"
  ) +
  coord_fixed()

print(p)
ggsave(file.path(dir.out, "axis_correlation_heatmap.png"),
       p, width = 5, height = 4.5, units = "in")
