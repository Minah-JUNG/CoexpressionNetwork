###############################################################
## 000. Module definition: PS / HEAT / ROS modules from subcluster IDs
###############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
})

## Module / subcluster pattern table
axis_tbl <- tribble(
  ~axis_code, ~axis_name,                                  ~pattern,
  "PS",       "PS-Energy Production & Maintenance",        "^C_1_4_|^C_1_1_|^C_6_1_1",
  "HEAT",     "Heat-Proteostasis & Cellular Protection",   "^C_6_2_|^C_6_3|^C_6_4|^C_1_2_2_",
  "ROS",      "ROS-Signal Integration & Defense",          "^C_1_3_1_1_1|^C_1_3_2_"
)

## Consistent palette used across all module figures.
axis_palette <- c(
  "PS"   = "#7FBF7F",
  "HEAT" = "#FB9A99",
  "ROS"  = "#edd072"
)

## Assign each row a module label based on its subcluster_id.

assign_module <- function(df, id_col = "subcluster_id") {
  pat_ps   <- axis_tbl$pattern[axis_tbl$axis_code == "PS"]
  pat_heat <- axis_tbl$pattern[axis_tbl$axis_code == "HEAT"]
  pat_ros  <- axis_tbl$pattern[axis_tbl$axis_code == "ROS"]

  df %>%
    mutate(
      axis = case_when(
        str_detect(.data[[id_col]], pat_ps)   ~ "PS",
        str_detect(.data[[id_col]], pat_heat) ~ "HEAT",
        str_detect(.data[[id_col]], pat_ros)  ~ "ROS",
        TRUE                                  ~ NA_character_
      )
    )
}
