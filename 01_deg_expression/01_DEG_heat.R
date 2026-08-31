# Heat stress DEGs: 25 vs 42 C at 1h and 1d timepoints (Nb leaf tissue)

##################################################

library(dplyr)
library(DESeq2)
library(matrixStats)   # rowMedians
library(tibble)

# === PATH ===
dir <- "/path/to"
# ============================================
setwd(dir)

##################################################

file_expr <- "./1_Data/NbHeat_ReadCount.txt"
data_expr <- read.table(file_expr, check.names = FALSE, header = T) 
data_expr <- data_expr %>%
  tibble::column_to_rownames(var = 'Geneid')

data <- data_expr[,!grepl("pc", colnames(data_expr))]

##------------------------------------------------
##################################################
# threshold setting

th.count <- 10
th.padj <- 0.05
th.fc <- 2

##################################################
## samples

samples <- data.frame(
  sample = colnames(data),
  Temp = ifelse(grepl("25", colnames(data)), "T25", "T42"),
  Time = ifelse(grepl("1D", colnames(data)), "D1", "H1")
)  %>%
  dplyr::mutate(condition = paste0(Temp, "_", Time))


samples$Temp <- factor(samples$Temp, levels = c("T25", "T42"))
samples$Time <- factor(samples$Time, levels = c("H1", "D1"))
samples$condition <- factor(samples$condition, levels = c("T25_H1", "T25_D1", "T42_H1", "T42_D1"))

##------------------------------------------------
## comparison group setting
##------------------------------------------------

samples
comparison <- samples$condition %>% unique
group_1 <- comparison[grepl("25", comparison)] %>% as.vector()
group_2 <- comparison[grepl("42", comparison)] %>% as.vector()


group_pairs <- data.frame(group_1, group_2)
group_pairs <- group_pairs[c(2,1),]

##------------------------------------------------

#group_pairs <- rbind(group_pairs.1, group_pairs.2)

##################################################

result_list <- list()

for(i in 1:nrow(group_pairs)){
  
  group_1 <- group_pairs[i,1]
  group_2 <- group_pairs[i,2]
  
  cat("Group 1 :", group_1, "\n")
  cat("Group 2 :", group_2, "\n")
  
  ##################################################
  # sample information
  
  samples_2pairs <- samples %>%
    filter(condition == group_1 | condition == group_2)
  samples_2pairs$condition <- factor(samples_2pairs$condition, levels=c(group_1, group_2))
  
  ##################################################
  # sample expression
  
  count <- data %>%
    dplyr::select(samples_2pairs$sample)
  
  # change
  count <- as.matrix(count) # for rowMedians
  count.ft <- count[rowMedians(count[,c(1:3)]) >= th.count | rowMedians(count[,c(4:6)]) >= th.count, ]
  
  ##################################################
  
  dds <- DESeqDataSetFromMatrix(countData = count.ft, colData = samples_2pairs, design = ~ condition)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  res <- res[order(res$padj),]
  
  ##################################################
  
  file.out <- paste0("DEGs_",group_1, "_", group_2,  "_", Sys.Date(),".txt")
  write.table(res,
              file = file.out, 
              row.names = T, 
              col.names = T, 
              sep = "\t")
  
  cat("===== DEGs result (individual):", file.out, "\n")
  
  ##################################################
  
  result_list[[i]] <- res
  names(result_list)[i] <- paste0(group_1, "_", group_2)
  
}

rm(file.out)

file.out <- paste0("./3_Results/Nicotiana_DEGs_Heat_", Sys.Date(), ".RData")
save(result_list,
     file = file.out)

cat("===== DEGs result (total):", file.out, "\n")

##################################################

th.fc <- 2

result_summary <- data.frame(
  Genes = sapply(result_list, nrow),
  Padj = sapply(result_list, function(x) subset(x, padj < th.padj)) %>% sapply(., nrow),
  FoldChange = sapply(result_list, function(x) subset(x, abs(log2FoldChange) > th.fc)) %>% sapply(., nrow),
  DEGs = sapply(result_list, function(x) subset(x, padj < th.padj & abs(log2FoldChange) > th.fc)) %>% sapply(., nrow),
  DEGs_UP = sapply(result_list, function(x) subset(x, padj < th.padj & log2FoldChange > th.fc)) %>% sapply(., nrow),
  DEGs_DW = sapply(result_list, function(x) subset(x, padj < th.padj & log2FoldChange < -(th.fc))) %>% sapply(., nrow)
)
print(result_summary)

rm(file.out)
file.out <- paste0("./3_Results/DEGs_HeatStress_summary_", Sys.Date(), ".txt")
write.table(result_summary, file.out, sep="\t")

cat("===== DEGs summary:", file.out, "\n")

##################################################

result_ft <- sapply(result_list, function(x) subset(x, padj < th.padj & abs(log2FoldChange) > th.fc))

##################################################
