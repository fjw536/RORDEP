library(DESeq2)
library(tidyverse)
library(ComplexHeatmap)
library(clusterProfiler)
library(ggrepel)

counts <- read.csv("counts.csv", row.names = 1)
metadata <- read.delim("metadata.txt", row.names = 1) %>%
  filter(sample_type == "liver")

dds <- DESeqDataSetFromMatrix(
  countData = counts[, rownames(metadata)],
  colData = metadata,
  design = ~ group
) %>%
  .[rowSums(counts(.) >= 5) >= ncol(.), ] %>%
  DESeq()

res <- lfcShrink(dds, coef = "group_r.RORDEP1_vs_Control", type = "ashr") %>%
  as.data.frame() %>%
  filter(abs(log2FoldChange) > 1 & padj < 0.05)

pca_data <- assay(vst(dds)) %>%
  t() %>%
  prcomp() %>%
  .$x %>%
  as.data.frame() %>%
  cbind(metadata)

ggplot(pca_data, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  labs(x = paste0("PC1 (", round(100 * summary(pca)$importance[2,1], "%)"))
       
       kegg <- enrichKEGG(
         gene = rownames(res),
         organism = "rno",
         pAdjustMethod = "BH",
         qvalueCutoff = 0.05
       )
       
       genes_of_interest <- c("Gene1", "Gene2")
       Heatmap(
         t(scale(t(assay(vst(dds)[genes_of_interest, ]))),
           name = "Z-score",
           show_row_names = TRUE,
           column_split = metadata$group
         )
         
         res %>%
           mutate(
             sig = case_when(
               padj < 0.05 & log2FoldChange > 1 ~ "Up",
               padj < 0.05 & log2FoldChange < -1 ~ "Down",
               TRUE ~ "NS"
             )
           ) %>%
           ggplot(aes(log2FoldChange, -log10(padj), color = sig)) +
           geom_point() +
           geom_text_repel(
             data = . %>% filter(sig != "NS"),
             aes(label = rownames(.)),
             max.overlaps = 20
           ) +
           geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
           geom_hline(yintercept = -log10(0.05), linetype = "dashed")
