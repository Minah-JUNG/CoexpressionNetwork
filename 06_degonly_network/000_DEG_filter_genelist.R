###############################################################
# 000. DEG filtering and gene-list extraction
###############################################################

library(dplyr)
library(tibble)
library(readr)

# === PATHS ===
# setwd("/path/to/working_dir")
# ===========================================

load("Nicotiana_DEGs_Heat_<DATE>_signicant.RData") # DEG from 02_deg_expression

## 1. 사용할 비교군
target_comps <- names(result_list)[1:2]
target_comps
# "T25_H1_T42_H1" "T25_D1_T42_D1"

## 2. 각 비교군에서 padj < 0.05 & abs(log2FC) > 2 유전자 추출
deg_list <- lapply(target_comps, function(comp) {
  as.data.frame(result_list[[comp]]) %>%
    rownames_to_column("Geneid") %>%
    filter(
      !is.na(padj),
      padj < 0.05,
      abs(log2FoldChange) > 2
    ) %>%
    mutate(comparison = comp)
})

names(deg_list) <- target_comps

## 3. 두 비교군 DEG union
deg_union <- unique(unlist(lapply(deg_list, function(x) x$Geneid)))

length(deg_union)

## 4. TPM 발현량 데이터 불러오기
data <- read.table(
  "Nb_ymkim_TPM.txt",
  header = TRUE,
  check.names = FALSE
)

## 5. DEG union 유전자만 TPM에서 선별
tpm_deg_total <- data %>%
  dplyr::filter(Geneid %in% deg_union)

## 저장: 전체 샘플 발현량
write.table(
  tpm_deg_total,
  file = "Nb_ymkim_TPM_DEG_TotalSample.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## 6. 샘플명에 25 또는 42가 들어간 heat 처리군 샘플만 선택
## H 샘플: 예) Nb-25-1H-A, Nb-42-1H-A 등
heat_cols <- grep("25.*H|42.*H|H.*25|H.*42", colnames(tpm_deg_total), value = TRUE)
heat_cols <- grep("-25-|-42-", colnames(tpm_deg_total), value = TRUE)

heat_cols
length(heat_cols)   # 12

tpm_deg_heat <- tpm_deg_total %>%
  dplyr::select(Geneid, all_of(heat_cols))

## 저장: heat 처리군 6개 샘플 발현량
write.table(
  tpm_deg_heat,
  file = "Nb_ymkim_TPM_DEG_HeatSample.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
