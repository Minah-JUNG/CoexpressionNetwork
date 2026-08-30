# Tissue-specific genes (one-vs-rest DESeq2 + CV-based background-stability filter)

library(dplyr)
library(DESeq2)
library(openxlsx)
library(matrixStats)   # rowMedians
library(tibble)

##################################################
# === PATH ===
dir <- "/path/to/Promoter_NbenthamianaTissue"
# ============================================
setwd(dir)

##################################################
file.data <- file.path(dir, "1_Data", "NbTissue_ReadCount.txt")

##################################################
# read.table(..., row.names = 1) already moves Geneid into rownames.
data <- read.table(file.data, header = TRUE, check.names = FALSE, row.names = 1)
colnames(data) <- stringr::str_replace_all(colnames(data), ".bam", "")

##################################################
# threshold setting
th.count <- 10
th.padj <- 0.05
th.fc <- 2

##################################################
# comparison group setting
samples <- data.frame(sample = colnames(data),
                      condition = c(rep("Flower",3),
                                    rep("Leaf",3),
                                    rep("Root",3),
                                    rep("Stem",3)))

group <- samples$condition %>% unique

##################################################
# main loop

for (target in unique(samples$condition)) {
  
  cat("Performing Tissue Specific Genes for:", target, "\n")
  
  samples_1vsrest <- samples %>%
    mutate(condition_new = ifelse(condition == target, target, "others"))
  samples_1vsrest$condition_new <- factor(samples_1vsrest$condition_new,
                                          levels = c("others", target))
  
  count_1vsrest <- data %>% dplyr::select(samples_1vsrest$sample)
  count_1vsrest <- as.matrix(count_1vsrest)
  
  target_idx <- which(samples_1vsrest$condition_new == target)
  others_idx <- which(samples_1vsrest$condition_new == "others")
  
  count.ft <- count_1vsrest[
    rowMedians(count_1vsrest[, target_idx]) >= th.count |
      rowMedians(count_1vsrest[, others_idx]) >= th.count, ]
  
  dds <- DESeqDataSetFromMatrix(
    countData = count.ft,
    colData = samples_1vsrest,
    design = ~ condition_new
  )
  dds <- DESeq(dds)
  
  # 정규화된 count
  norm_counts <- counts(dds, normalized=TRUE)
  
  # DEG 결과 (target vs rest)
  res <- results(dds)
  res <- res[order(res$padj), ]
  res_df <- as.data.frame(res)
  
  # 정규화된 count
  norm_counts_df <- as.data.frame(norm_counts)
  norm_counts_df$gene <- rownames(norm_counts_df)
  
  # DEG 결과
  res_df$gene <- rownames(res_df)
  
  # (output1) rowname 기준으로 merge
  output1 <- merge(norm_counts_df,
                   res_df[, c("gene","log2FoldChange","padj")],
                   by = "gene", all.x = TRUE)
  rownames(output1) <- output1$gene
  output1$gene <- NULL
  
  output2 <- subset(output1, padj < th.padj & log2FoldChange > th.fc)
  output2 <- output2[order(output2$padj), ]
  
  #----------------------------------------------
  # others 그룹 내 분산 기준 filtering
  #----------------------------------------------
  
  others_samples <- samples %>% filter(condition != target)
  count_others <- data %>% dplyr::select(others_samples$sample) %>% as.matrix()
  
  # 유전자별 평균, 분산 계산
  gene_mean <- rowMeans(count_others)
  gene_var  <- apply(count_others, 1, var)
  
  # CV 
  gene_cv <- sqrt(gene_var) / (gene_mean + 1e-8)
  
  # 분산/변동계수 기준 필터링
  th.cv <- 0.5
  bg_stable_genes <- names(gene_cv[gene_cv < th.cv])
  cat("Background stable genes (variance filter):", length(bg_stable_genes), "\n")
  
  # 최종 tissue-specific 후보군
  final_genes <- intersect(rownames(output2), bg_stable_genes)
  output3 <- output1[final_genes, ]
  
  write.table(output3,
              paste0("./3_Results/Nb_TissueSpecific_", target, "_up_VarianceFiltered_DESeq2.csv"),
              sep = ",", quote = FALSE, row.names = TRUE)
  
}


