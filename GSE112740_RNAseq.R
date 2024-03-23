
setwd("E:/Master_3rd_year/Projects/Guanchuan_RNAseq")

library(DESeq2)

#### 导入RSEM.results表格
f.ls <- list.files("GSE112740_RAW", pattern = "txt$", full.names = T)

for (i in 1:length(f.ls)) {
  rsem_data <- read.table(list.files(f.ls, full.names = T)[i], header = TRUE, row.names = 1)
  colnames(rsem_data)
  # rsem_data <- rsem_data[,c("transcript_id.s.", "pme_TPM")]
  rsem_data <- rsem_data[,c("transcript_id.s.", "expected_count")]
  
  rsem_data$gene <- lapply(strsplit(rownames(rsem_data), "_"), function(x){x[[2]]}) %>% unlist
  # 简单粗暴
  rsem_data <- rsem_data[!duplicated(rsem_data$gene),]
  rownames(rsem_data) <- rsem_data$gene
  rsem_data <- rsem_data[,-1]
  colnames(rsem_data)[1] <- (list.files(f.ls)[i] %>% strsplit(split = "\\.RSEM"))[[1]][1]
  head(rsem_data)
  
  if (i==1) {
    rsem_mtx <- rsem_data
  } else {
    rsem_mtx <- left_join(rsem_mtx, rsem_data)
  }
}
rownames(rsem_mtx) <- rsem_mtx$gene
rsem_mtx <- rsem_mtx[,-2]
head(rsem_mtx)


#### 创建样本信息数据框
sample_info <- data.frame(
  sample = colnames(rsem_mtx),
  condition = c("HighFat", "HighFat", "HighFat", "LowFat", "LowFat", "LowFat")
)
sample_info


#### 创建DESeq2对象
rsem_mtx.int <- round(rsem_mtx)
dds <- DESeqDataSetFromMatrix(countData = rsem_mtx.int, colData = sample_info, design = ~ condition)
#### 运行差异表达分析
dds <- DESeq(dds)
#### 提取差异基因
res <- results(dds)
write.csv(res, "Results2.csv")






