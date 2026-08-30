# Gene-pair correlation matrix from TPM data.
# Pipeline: Mock excluded -> log2(TPM+1) transform -> Pearson correlation.

rm(list = ls())

library(dplyr)
library(parallel)
library(Matrix)
library(tibble)
library(arrow)

# === PATH ===
# setwd("/path/to/correlation_input_dir")
# ===========================================

# -------------------------------------------
# Load TPM data
# -------------------------------------------
file <- "Nb_ymkim_TPM.txt" # merge NbHeat and NbTissue
data.nor <- read.table(file, check.names = FALSE, header = TRUE)
data.nor <- data.nor %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("Geneid")

# 전치 (행: 샘플, 열: 유전자)
data_ft <- data.nor[, colSums(data.nor) != 0]  # all-zero 제거
data_ft <- as.matrix(data_ft) %>% t

cat("Samples:", nrow(data_ft), " Genes(before filter):", ncol(data_ft), "\n")

# -------------------------------------------
# 필터링 조건 적용 (3Q=0 & rowMax <=10)
# -------------------------------------------
q3_vals   <- apply(data_ft, 2, quantile, probs=0.75, na.rm=TRUE)
row_max   <- apply(data_ft, 2, max, na.rm=TRUE)

exclude_idx <- (q3_vals == 0 & row_max <= 10)
include_idx <- !exclude_idx

cat("Excluded genes:", sum(exclude_idx), "\n")
cat("Remaining genes:", sum(include_idx), "\n")

# 제외된 TPM 데이터 저장
excluded_tpm <- data_ft[, exclude_idx, drop=FALSE]
saveRDS(excluded_tpm, "Excluded_genes_TPM.rds")
write.csv(excluded_tpm, "Excluded_genes_TPM.csv")

# 필터링된 데이터만 사용
data_ft <- data_ft[, include_idx, drop=FALSE]

# 변동 없는 유전자 제거
data_ft <- data_ft[, apply(data_ft, 2, sd) > 0]

cat("Final samples:", nrow(data_ft), " Genes(after filter):", ncol(data_ft), "\n")

# -------------------------------------------
# log2 transform (TPM → log2(TPM+1))
# -------------------------------------------
data_ft <- log2(data_ft + 1)

# -------------------------------------------
# 전체 상관계수 계산
# -------------------------------------------
cor_mat <- cor(data_ft, method="pearson", use="pairwise.complete.obs")

cat("===== Done: correlation calculated ===== \n")

# -------------------------------------------
# p-value 계산 함수
# -------------------------------------------
cor2p <- function(r, n) {
  tval <- r * sqrt((n-2) / (1-r^2))
  2 * pt(-abs(tval), df = n-2)
}

n <- nrow(data_ft)
p_mat <- cor2p(cor_mat, n)

cat("===== Done: p-value calculated ===== \n")

# -------------------------------------------
# 결과 저장 (sparse matrix 권장)
# -------------------------------------------
cor_sparse <- Matrix(cor_mat, sparse=TRUE)
p_sparse   <- Matrix(p_mat, sparse=TRUE)

saveRDS(cor_sparse, "Correlation_matrix_sparse.rds")
saveRDS(p_sparse, "Pvalue_matrix_sparse.rds")

cat("Done! Results saved as sparse matrix.\n")

# -------------------------------------------
# 필터링 함수 정의
# -------------------------------------------
get_filtered_pairs <- function(cor_mat, p_mat, thr_r, thr_p, chunk_size = 1e6) {
  upper_idx <- which(upper.tri(cor_mat), arr.ind = TRUE)
  total <- nrow(upper_idx)
  
  results <- list()
  for (i in seq(1, total, by = chunk_size)) {
    j <- min(i + chunk_size - 1, total)
    
    sel <- abs(cor_mat[upper_idx[i:j, , drop=FALSE]]) > thr_r &
           p_mat[upper_idx[i:j, , drop=FALSE]] < thr_p
    
    if (any(sel)) {
      results[[length(results) + 1]] <- data.frame(
        Gene1 = colnames(cor_mat)[upper_idx[i:j, 1]][sel],
        Gene2 = colnames(cor_mat)[upper_idx[i:j, 2]][sel],
        Correlation = cor_mat[upper_idx[i:j, , drop=FALSE]][sel],
        Pvalue = p_mat[upper_idx[i:j, , drop=FALSE]][sel],
        stringsAsFactors = FALSE
      )
    }
  }
  
  dplyr::bind_rows(results)
}

# -------------------------------------------
# Supplementary 결과 출력
# -------------------------------------------
thresholds <- c(0.70, 0.75, 0.80, 0.85, 0.90, 0.95)
threshold_p <- 0.05

for (thr_r in thresholds) {
  cat("\n--- Threshold |r| >", thr_r, " & p <", threshold_p, " ---\n")
  t0 <- Sys.time()
  filtered_pairs <- get_filtered_pairs(cor_mat, p_mat, thr_r, threshold_p)
  t1 <- Sys.time()
  
  cat("Elapsed time:", round(as.numeric(difftime(t1, t0, units="secs")),2), "sec\n")
  cat("Related gene pairs:", nrow(filtered_pairs), "\n")
  cat("Unique genes:", length(unique(c(filtered_pairs$Gene1, filtered_pairs$Gene2))), "\n")
}

# -------------------------------------------
# 최종 저장 (r=0.7, p<0.05 기준)
# -------------------------------------------
threshold_r <- 0.7
filtered_pairs <- get_filtered_pairs(cor_mat, p_mat, threshold_r, threshold_p)

saveRDS(filtered_pairs, "Filtered_correlation_pairs.rds")
write.csv(filtered_pairs, "Filtered_correlation_pairs.csv", row.names = FALSE)

# parquet 저장은 맨 뒤 (매우 느림)
write_parquet(filtered_pairs, "Filtered_correlation_pairs.parquet")

cat("Results saved: RDS, CSV, and Parquet formats.\n")
