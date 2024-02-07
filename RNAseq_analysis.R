library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(data.table)
library(ggThemeAssist)
library(dplyr)
library(ggpubr)
library(UpSetR)
library(clusterProfiler)
library(ReactomePA)
library(biomaRt)
library(enrichplot)
library(DESeq2)
library(pcaExplorer)
library(apeglm)
library(ggplot2)
library(ggfortify)
library(patchwork)
library(conflicted)
library(writexl)
library(ashr)
library(openxlsx)
library(mixOmics)
library(RColorBrewer)
library(ggvenn)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(scales)
library(ggridges)
library(GOSemSim)
library(readxl)
library(DOSE)
library(forcats)
library(msigdbr)
library(KEGGREST)
library(ggplotify)
library(patchwork)
library(ggpubr)
library(scales)
library(edgeR)
library(pheatmap)
library(viridis)
library(ggrepel)
library(KEGGREST)
library(factoextra)

# setwd("/Users/nhx573/Library/CloudStorage/OneDrive-UniversityofCopenhagen/work_in_2020-2023/RUCILP/bulk_RNAseq_start20230126/bulkRNAseq_20231228")
genecounts_all <- read.delim("counts_20231228_geneID.csv", sep = ",", header = T, row.names=1) 
annotation <- read.delim("Annotation.txt", header = T, row.names = 1)
out_dir <- "liver"
tissue_type <- "liver"
metadata_all <- read.delim("metadata_all20231228.txt", header = T, row.names = 1)
metadata <- metadata_all %>% clusterProfiler::filter(sample_type == tissue_type)
genecounts <- as.data.frame(apply(genecounts_all[, rownames(metadata)], 2, as.integer))
rownames(genecounts) <- rownames(genecounts_all[, rownames(metadata)])
metadata$group <- factor(metadata$group, levels = c("Control", "r-RORDEP1"))
dds <- DESeqDataSetFromMatrix(countData = genecounts,
                              colData = metadata,
                              design = ~ group)
genecounts_keep <- rowSums(counts(dds) >= 5) >= nrow(metadata)
dds_filtered <- dds[genecounts_keep,]
dds_filtered <- DESeq(dds_filtered)
sizefactors <- sizeFactors(dds_filtered)
library_sizes <- colSums(counts(dds_filtered))
dds_filtered_normalizedcounts <- as.data.frame(counts(dds_filtered, normalized=TRUE), row.names = names(dds_filtered))

list_res <- list()
list_res_filter <- list()
list_lfcShrink <- list()
list_lfcShrink_filter <- list()

for(i in 1:(length(all_groups) - 1)) {
  for(j in (i + 1):length(all_groups)) {
    group1 <- all_groups[i]
    group2 <- all_groups[j]
    contrast <- c("group", group2, group1)
    res_custom = results(dds_filtered, contrast=contrast)
    res_custom$annotation = annotation$gene_name[match(rownames(res_custom), rownames(annotation))]
    filtered_res_custom = subset(res_custom, padj < 0.05 & abs(log2FoldChange) > log2(2))
    res_custom_shrink = lfcShrink(dds_filtered, res=res_custom, type="ashr")
    res_custom_shrink$annotation = annotation$gene_name[match(rownames(res_custom_shrink), rownames(annotation))]
    filtered_res_custom_shrink = subset(res_custom_shrink, padj < 0.05 & abs(log2FoldChange) > log2(2))
    list_res[[paste0(group2, "_vs_", group1)]] <- res_custom
    list_res_filter[[paste0(group2, "_vs_", group1)]] <- filtered_res_custom
    list_lfcShrink[[paste0(group2, "_vs_", group1)]] <- res_custom_shrink
    list_lfcShrink_filter[[paste0(group2, "_vs_", group1)]] <- filtered_res_custom_shrink
  }
}

group_colors <- setNames(palette_GR_BL, c("Control", "r-RORDEP1"))
output_dds_filtered_vst <- as.data.frame(assay(vst(dds_filtered, blind = TRUE)))
pca_result_vst <- prcomp(t(output_dds_filtered_vst), center = TRUE, scale. = TRUE)
percentage_variance_explained_vst <- summary(pca_result_vst)$importance["Proportion of Variance", ] * 100
pc_scores_vst <- as.data.frame(pca_result_vst$x)
pca_data_vst <- cbind(metadata, pc_scores_vst)
sampleID <- rownames(pca_data_vst)
group_order <- c("Control","r-RORDEP1")
pca_data_vst$group <- factor(pca_data_vst$group, levels = group_order)
PCA_vst_plot <- ggplot(pca_data_vst, aes(x = PC1, y = PC2, colour = group, shape = group)) +
  geom_point(size = 10, alpha = 0.7, position = position_jitter(width = 2, height = 2))+
  scale_color_manual(values = group_colors) +
  coord_fixed(xlim = c(-100, 100), ylim = c(-100, 100)) +
  scale_x_continuous(limits=c(-100, 100)) +  
  scale_y_continuous(limits=c(-100, 100)) + 
  labs(
    x = paste0("PC1 (", round(percentage_variance_explained_vst[1], 0), "%)"),
    y = paste0("PC2 (", round(percentage_variance_explained_vst[2], 0), "%)")) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text( size = 26),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1, 1, 1, 1, unit = "cm"),
    axis.title = element_text(margin = margin(r = 12), size = 29),
    axis.text = element_text(size = 26)
  ) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted")

perform_all_enrichments <- function(df, out_dir, tissue_type, direction) {
  input_gene = row.names(df)
  KEGG <- enrichKEGG(
    gene = input_gene,
    organism = 'rno',
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 8,
    maxGSSize = 500
  )
  KEGG <- setReadable(KEGG,
                      OrgDb = org.Rn.eg.db,
                      keyType = "ENTREZID")
  
  enrichments <- list(KEGG = KEGG)
  wb <- createWorkbook()
  for (name in names(enrichments)) {
    addWorksheet(wb, name)
    writeData(wb, name, enrichments[[name]])
  }
  saveWorkbook(wb, paste0(out_dir, '/', paste0(direction, "_Enrichment_Analysis_", tissue_type, ".xlsx")), overwrite = TRUE)
  return(enrichments)
}

original_DEGs <- as.data.frame(list_lfcShrink_filter$"r-RORDEP1_vs_Control")
upregulated_DEGs <- subset(original_DEGs, log2FoldChange > 1)
downregulated_DEGs <- subset(original_DEGs, log2FoldChange < -1)

direction <- c("Original", "Upregulated", "Downregulated")
deg_lists <- list(original_DEGs, upregulated_DEGs, downregulated_DEGs)

enrich_list <- list()
for (i in 1:length(direction)) {
  enrich_list[[direction[i]]] <- perform_all_enrichments(deg_lists[[i]], out_dir, tissue_type, direction[i])
}



focusofsubtype <- c("Original", "subcategory", "category")
enrich_category <- enrich_list[[focusofsubtype[1]]][["KEGG"]]@result
enrich_category <- enrich_category[enrich_category$p.adjust < 0.05, ]
category_data <- as.data.frame(table(enrich_category[, c( focusofsubtype[3], focusofsubtype[2])]))
names(category_data) <- c( "category",  "subcategory", 'count')
category_data <- category_data[category_data$count>0, ]

category_data <- category_data %>%
  transform(percentage = round((count / sum(count)) * 100, 1)) %>%
  arrange(desc(count)) %>%
  mutate(label = paste(category, "\n", percentage, "%"))

category_data <- category_data %>%
  mutate(category = factor(category, levels = unique(category[order(-count)])))

category_data_2_modified <- category_data_2 %>%
  mutate(Subcategory = ifelse(count <= 1, "Others", Subcategory)) %>%
  group_by(Subcategory) %>%
  summarise(count = sum(count), .groups = 'drop')

category_data_2_modified <- category_data_2_modified %>%
  transform(percentage = round((count / sum(count)) * 100, 1)) %>%
  arrange(desc(count)) %>%
  mutate(label = ifelse(row_number() <= 6, paste(Subcategory, "\n", percentage, "%"), ""))

category_data_2_modified$label[category_data_2_modified$Subcategory == "Others"] <- ""

enrich_donut <- ggplot(category_data_2_modified, aes(x = "", y = percentage, fill = Subcategory)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  theme_void() +
  theme(legend.position = "left") +
  annotate("text", x = 0, y = 0, label = paste0("85 enriched pathways", "\n", "by KEGG"), size = 4.5) +
  geom_text(aes(label = str_wrap(label, width = 20)), position = position_stack(vjust = 0.5), size = 4)

dev.off()

mygenes2_path <- paste0(out_dir, '/mygenes_', tissue_type, '_2.xlsx')
sheet_names <- excel_sheets(mygenes2_path)
process_sheet <- function(sheet_name, file_path, annotation) {
  df <- read_excel(file_path, sheet = sheet_name)
  gene_id <- row.names(annotation)[match(df$gene_symbol, annotation$gene_name)]
  df <- transform(df, gene_ID = gene_id)
  
  df <- df[!is.na(df$gene_ID), ]
  return(df)
}

mygenes2 <- setNames(lapply(sheet_names, process_sheet, file_path = mygenes2_path, annotation = annotation), sheet_names)
mygene2_symbols <- unlist(lapply(mygenes2, function(df) df$gene_symbol))
unique_mygene2_symbols <- unique(mygene2_symbols)

mygenes2_df <- na.omit(allin_df[allin_df$annotation %in% unique_mygene2_symbols & allin_df$padj < 0.05, ])
mygenes2_df <- mygenes2_df[, c(1:12)]
mygenes2_df_z_score <- as.data.frame(t(scale(t(mygenes2_df))))
mygenes2_df_z_score <- `rownames<-`(mygenes2_df_z_score, annotation$gene_name[match(row.names(mygenes2_df_z_score), rownames(annotation))])
unique_mygene2_symbols <- factor(unique_mygene2_symbols, levels = unique_mygene2_symbols)
mygenes2_df_z_score <- na.omit(mygenes2_df_z_score[match(levels(unique_mygene2_symbols), row.names(mygenes2_df_z_score)), ])

gene_category_color <- c("#008080","#016064", "#3F704D", "#c35b3b", "#cc7154","orange")
gene_category <- c("Gluconeogenesis","Glycogenolysis", "Lipogenesis","Glycogenesis","Glycolysis" ,"Insulin-PI3K-AKT signaling")
ha = rowAnnotation(foo = anno_empty(border = F, width = unit(1.5, "cm")))

mygenes2_heatmap2 <- ComplexHeatmap::Heatmap(
  mygenes2_df_z_score, 
  width = unit(5, "cm"),
  name = "mygenes2_heatmap",
  col = heatmap_color,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  show_row_names = T,
  show_column_names = F,
  show_row_dend = F,
  show_column_dend = T, 
  border = F,  
  top_annotation = group_annotation,
  right_annotation = ha, 
  column_title=NULL,
  row_title_rot = 0,
  row_names_side = "left",
  row_names_rot = 0,
  row_title= NULL,
  row_names_gp = gpar(fontsize = 8, fontface = "italic"),
  row_title_gp = gpar(fontsize = 8),
  show_heatmap_legend =F,
  row_split = gene_annotation_df$gene_category,
  column_split = group_annotation_df$group,
  column_gap = unit(0.1, "mm"),
  row_gap = unit(1, "mm"),
  column_dend_height = unit(1, "cm")
)

png(paste0(out_dir, '/', tissue_type, '_interestedgenes_beneficial_complexHM3.png'), width = 2350, height = 3800, res = 600)
draw(mygenes2_heatmap2, merge_legend = TRUE, padding = unit(c(5, 5, 5, 5), "mm"))
for(i in 1:6) {
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, width = unit(0.5, "mm"), gp = gpar(fill = "black", col = NA), just = "left")
    grid.text(paste(strwrap(gene_category[i], width = 2), collapse = "\n"), x = unit(3, "mm"), just = "left", gp = gpar(fontfamily = "Arial", fontsize = 8))
  })
}
dev.off()

threshold_padj <- 0.05
lfc_threshold <- log2(2)
color_volcano <- c("Downregulated DEGs" = "#f09e83", "Upregulated DEGs" = "#346fc4", "Non-significant" = "#b5b5b5")
volcano <- original_DEGs <- as.data.frame(list_lfcShrink$"r-RORDEP1_vs_Control")
volcano$sig <- with(volcano, ifelse(log2FoldChange < -lfc_threshold & padj < threshold_padj, "Downregulated DEGs",
                                    ifelse(log2FoldChange > lfc_threshold & padj < threshold_padj, "Upregulated DEGs", "Non-significant")))
volcano$highlight <- ifelse(volcano$annotation %in% row.names(mygenes2_df_z_score), "Highlight", "NoHighlight")
highlighted_data <- subset(volcano, highlight == "Highlight")
volcano <- na.omit(volcano)
volcano$sig <- factor(volcano$sig, levels = c("Upregulated DEGs", "Downregulated DEGs", "Non-significant"))
scale_color_manual(values = color_volcano, 
                   name = NULL,
                   breaks = c("Upregulated", "Downregulated"))
volcano.plot2 <- ggplot(volcano, aes(x=log2FoldChange, y=-log10(padj), color=sig)) +
  geom_point(alpha=0.4) +
  geom_text_repel(size = 5, fontface = "italic", data=highlighted_data,
                  aes(x=log2FoldChange, y=-log10(padj), label=annotation),
                  segment.color = 'black', color = 'black') +
  theme_minimal() +
  labs(x = "Log2 fold change", y = "-Log10 adjusted p value") +
  scale_x_continuous(breaks=seq(-5, 5, by=2), limits=c(-5.1, 5.1)) +
  scale_y_continuous(breaks=seq(0, 50, by=15), limits=c(0, 48)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        axis.text = element_text(size = 26), 
        axis.title = element_text(size = 28),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size =20)) +
  scale_color_manual(values = color_volcano, 
                     name = NULL,
                     breaks = c("Upregulated DEGs", "Downregulated DEGs")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  geom_hline(yintercept = -log10(threshold_padj), color="gray", linetype="dotted") +
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), color="gray", linetype="dotted") +
  annotate("text", x = -5.1, y = -log10(0.001), label = paste("padj=", threshold_padj), hjust = 0, vjust = 1.1, size = 5, fontface = "italic", color = "black")
print(volcano.plot2)
dev.off()
ggsave(paste0(out_dir, '/', tissue_type, '_volcano2.pdf'), plot = volcano.plot2, width = 8, height = 8, dpi = 300)
