#load libraries
library(GOplot)
library(reshape2)
library(pheatmap)
library(ggrepel)
library(limma)
library(robustbase)  
library(DESeq2)
library(heatmap3)
library(RColorBrewer)
library(clusterProfiler)
library(stringr)
library(EnhancedVolcano)
library(org.Mm.eg.db)
library(ggplot2)
library(ComplexHeatmap)
library(tidyr)
library(dplyr)
library(grid)
library(gtable)
library(gridGraphics)
library(ggrepel)
library(openxlsx)
library(pathview)
library(enrichplot)
library(circlize)
library(tidyverse)



getwd()
setwd("Gabrg3Bulk")
#load Data
data <- read.csv("gene_expression_rawcounts.csv", header=T, row.names="gene_symbol")
#data <- data[,1:12]
meta <- data.frame(
  sample = colnames(data),
  condition = rep(c("WT","KO"),each=6),
  region = rep(rep(c("ctx","hp"),each=3),2)
)
meta <- meta %>% mutate(
  condition = factor(condition, levels= c("WT","KO")),
  region = factor(region, levels= c("ctx","hp"))
)

count_data <- apply(data, 2, function(x) {
  round(as.numeric(str_extract(x, "\\d+\\.?\\d*")))
})
rownames(count_data) <- rownames(data)


design_formula <- ~ region + condition + region:condition
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = meta, design = design_formula)
dds <- DESeq(dds)
resultsNames(dds)

res_inter <- results(dds, name="regionhp.conditionKO") 

res_ctx <- results(dds, contrast=c("condition", "KO", "WT"))
res_hp <- results(dds, contrast=list(c("condition_KO_vs_WT", "regionhp.conditionKO")))

res_ctx_sig <- subset(res_ctx, padj < 0.05 & abs(log2FoldChange) > 1)
res_hp_sig <- subset(res_hp, padj < 0.05 & abs(log2FoldChange) > 1)

res_ctx_sig1.5 <- subset(res_ctx, padj < 0.05 & abs(log2FoldChange) > 0.585)
res_hp_sig1.5 <- subset(res_hp, padj < 0.05 & abs(log2FoldChange) > 0.585)

res_c <- as.data.frame(res_ctx)
res_h <- as.data.frame(res_hp)
shown_genes_ctx <- c("Esr1", "Ptgs2", "Epha4", "Neu4", "Mc4r", "Draxin", "Fos", "Postn", "Gper1", "Nrn1", "Gad2", "Rbfox3", "Lzts1", "Gabra1", "Dlx2", "Arx", "Dlx1", "Lrtm2")
shown_genes_hp <- c("Esr1", "Dcn", "Nlrp3", "Serpinf1", "Penk", "Serpine1", "Gria1", "Itgam", "Aqp4", "Ugt1a7c", "Mertk", "Slc7a11", "Trem2", "Spp1", "Dock2", "Ptk2b", "Lpl", "Pirb", "Tyrobp", "Cd36")

gene <- c("Esr1", "Pdgfra", "Col1a1", "Col3a1", "Col4a1", "Fn1", "Thbs1", "Itgav", "Il1rap", "Dusp6", "Spry1", "Spry4", "Mbp", "Gfap","Nrxn3","Nrxn1","Nlgn1","Gphn")

shown_ctx <- subset(res_ctx, rownames(res_ctx) %in% union(shown_genes_ctx,gene))
shown_hp <- subset(res_hp, rownames(res_ctx) %in% union(shown_genes_hp,gene))

genes_c <- as.data.frame(res_ctx_sig[rownames(res_ctx_sig) %in% shown_genes_ctx,])
genes_h <- as.data.frame(res_hp_sig[rownames(res_hp_sig) %in% shown_genes_hp,])

cols_c <- ifelse(res_c$padj < 0.05 & res_c$log2FoldChange > 1, 'red', 
                 ifelse(res_c$padj < 0.05 & res_c$log2FoldChange < -1, 'blue',
                        'grey80')) 
cols_c[is.na(cols_c)] <- 'grey80'
names(cols_c)[cols_c == 'red'] <- 'Upregulated (KO > WT)'
names(cols_c)[cols_c == 'blue'] <- 'Downregulated (WT > KO)'
names(cols_c)[cols_c == 'grey80'] <- 'Non-Significant'

cols_h <- ifelse(res_h$padj < 0.05 & res_h$log2FoldChange > 1, 'red', 
                 ifelse(res_h$padj < 0.05 & res_h$log2FoldChange < -1, 'blue',
                        'grey80')) 
cols_h[is.na(cols_h)] <- 'grey80'
names(cols_h)[cols_h == 'red'] <- 'Upregulated (KO > WT)'
names(cols_h)[cols_h == 'blue'] <- 'Downregulated (WT > KO)'
names(cols_h)[cols_h == 'grey80'] <- 'Non-Significant'

## volcano plot
# cortex
EnhancedVolcano(res_c,
                lab = rownames(res_c),
                x = "log2FoldChange",
                y = "padj",
                colCustom = cols_c,
                colAlpha = 0.8,
                pCutoff = 0.05, 
                FCcutoff = 1,
                selectLab = genes_c,
                title = "Cortex (KO v.s. WT)",
                subtitle= "Differentially Expressed Genes")
# hip
EnhancedVolcano(res_h,
                lab = rownames(res_h),
                x = "log2FoldChange",
                y = "padj",
                colCustom = cols_h,
                colAlpha = 0.8,
                pCutoff = 0.05, 
                FCcutoff = 1,
                selectLab = rownames(res_hp_sig),
                title = "Hippocampus (KO v.s. WT)",
                subtitle= "Differentially Expressed Genes")


# Venn plot

frame_ctx <- as.data.frame(res_ctx_sig)
frame_hp <- as.data.frame(res_hp_sig)

frame_ctx <- frame_ctx %>% mutate(gene = rownames(frame_ctx))
frame_ctx <- select(frame_ctx, c(gene, log2FoldChange))
frame_hp <- frame_hp %>% mutate(gene = rownames(frame_hp))
frame_hp <- select(frame_hp, c(gene, log2FoldChange))

venn <- GOVenn(frame_ctx,frame_hp, label = c('KO in Cortex', 'KO in Hip'))
venn$plot_env$table


# scatter plot
res_df <- data.frame(
  gene = rownames(res_ctx),
  lfc_ctx = res_ctx$log2FoldChange,
  lfc_hp = res_hp$log2FoldChange,
  padj_ctx = res_ctx$padj,
  padj_hp = res_hp$padj,
  padj_int = res_inter$padj
)

### IMPORTANT: whether to consider the padj?
res_df <- res_df %>% mutate(
  Group = case_when(
    ((lfc_hp > 0 & lfc_ctx < 0) & (padj_ctx < 0.05 & padj_hp < 0.05)) ~ "Hp_up_Ctx_down",
    ((lfc_hp < 0 & lfc_ctx > 0) & (padj_ctx < 0.05 & padj_hp < 0.05)) ~ "Hp_down_Ctx_up",
    ((lfc_hp > 0 & lfc_ctx > 0 ) & (padj_ctx < 0.05 & padj_hp < 0.05)) ~ "Both_up",
    ((lfc_hp < 0 & lfc_ctx < 0 ) & (padj_ctx < 0.05 & padj_hp < 0.05)) ~ "Both_down",
    TRUE ~ "Non_Significant"
  )
)
Group_colors <- c(
  'Hp_up_Ctx_down' = '#984EA3',
  'Hp_down_Ctx_up' = '#4DAF4A',
  'Both_up' = '#E41A1C',
  'Both_down' = '#6495ED',
  'Non_Significant' = 'grey80'
)

# shared up/down-regulated genes
res_plot <- res_df %>%
  filter(gene %in% union(rownames(res_ctx_sig), rownames(res_hp_sig))) %>%
  mutate(dist = sqrt(lfc_ctx^2 + lfc_hp^2)) 

top_genes <- res_plot %>%
  slice_max(dist, n = 10) 

# top_hits <- res_df %>%
#   filter(padj_int < 0.05) %>% filter(Group %in% c("Hp_up_Ctx_down","Hp_down_Ctx_up")) %>%
#   group_by(Group) %>%
#   slice_min(padj_int, n = 5) 

# opposite genes
top_ctx <- res_df %>%
  filter(Group %in% c("Hp_up_Ctx_down","Hp_down_Ctx_up")) %>%
  mutate(dist = lfc_ctx - lfc_hp) %>%
  group_by(Group) %>%
  slice_max(dist, n = 5) 

top_hp <- res_df %>%
  filter(Group %in% c("Hp_up_Ctx_down","Hp_down_Ctx_up")) %>%
  mutate(dist = lfc_hp - lfc_ctx) %>%
  group_by(Group) %>%
  slice_max(dist, n = 5) 

# genes shown
genes <- union(top_genes$gene,union(top_ctx$gene, top_hp$gene))

ggplot(res_df, aes(x = lfc_ctx, y = lfc_hp)) +
  geom_point(data = res_df, aes(color = Group), size = 1.5, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_text_repel(data = filter(res_df, gene %in% genes),
                  aes(label = gene),
                  box.padding = 0.5, max.overlaps = Inf) +
  theme_classic() +
  scale_color_manual(values = Group_colors) +
  labs(title = "Differential Response: Ctx vs Hp",
       x = "Log2 Fold Change (mPFC)", y = "Log2 Fold Change (Hip)")

# GO bubble

gene_list <- list(
  Shared_Up = res_df$gene[res_df$Group == "Both_up"],
  Shared_Down = res_df$gene[res_df$Group == "Both_down"],
  Hp_up_Ctx_down = res_df$gene[res_df$Group == "Hp_up_Ctx_down"],
  Hp_down_Ctx_up = res_df$gene[res_df$Group == "Hp_down_Ctx_up"]
)

ck <- compareCluster(geneCluster = gene_list, fun = "enrichGO", 
                     OrgDb = org.Mm.eg.db, ont = "BP", keyType = "SYMBOL")
dotplot(ck) 

go_res <- as.data.frame(ck@compareClusterResult)

write.xlsx(go_res,"go_res.xlsx")

# for shared_up/down

shared_hp <- mapIds(org.Mm.eg.db, keys = res_df$gene[res_df$Group == "Both_up"], column = "ENTREZID", keytype = "SYMBOL")

go_shared_hp <- enrichGO(shared_hp, 
                         OrgDb = org.Mm.eg.db, 
                         ont = "BP",
                         pAdjustMethod = "BH",readable = TRUE,
                         pvalueCutoff = 0.05)

shared_down <- mapIds(org.Mm.eg.db, keys = res_df$gene[res_df$Group == "Both_down"], column = "ENTREZID", keytype = "SYMBOL")

go_shared_down <- enrichGO(shared_down, 
                           OrgDb = org.Mm.eg.db, 
                           ont = "BP",
                           pAdjustMethod = "BH",readable = TRUE,
                           pvalueCutoff = 0.05)

go_shared_down_filtered <- subset(go_shared_down@result, p.adjust < 0.05)
go_shared_hp_filtered <- subset(go_shared_hp@result, p.adjust < 0.05)
saveRDS(go_shared_hp,file = "go_shared_hp.rds") 
saveRDS(go_shared_down,file = "go_shared_down.rds") 
write.xlsx(go_shared_down_filtered[order(go_shared_down_filtered$Count, decreasing = T),],"go_shared_down.xlsx")
write.xlsx(go_shared_hp_filtered[order(go_shared_hp_filtered$Count, decreasing = T),],"go_shared_hp.xlsx")

# for res_inter
res_inter_sig <- subset(res_inter, padj < 0.05 & abs(log2FoldChange) > 0.585) 
res_inter_up <- subset(res_inter_sig, log2FoldChange > 0)
res_inter_down <- subset(res_inter_sig, log2FoldChange < 0)

inter_up <- mapIds(org.Mm.eg.db, keys=rownames(res_inter_up), column="ENTREZID", keytype="SYMBOL")
go_inter_up <- enrichGO(
  gene = inter_up,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH"
)
inter_down <- mapIds(org.Mm.eg.db, keys=rownames(res_inter_down), column="ENTREZID", keytype="SYMBOL")
go_inter_down <- enrichGO(
  gene = inter_down,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH"
)

go_inter_down_filtered <- subset(go_inter_down@result, p.adjust < 0.05) 
go_inter_up_filtered <- subset(go_inter_up@result, p.adjust < 0.05) 

write.xlsx(go_inter_down_filtered[order(go_inter_down_filtered$Count,decreasing = T),],"go_inter_down.xlsx")
write.xlsx(go_inter_up_filtered[order(go_inter_up_filtered$Count,decreasing = T),],"go_inter_up.xlsx")

go_combined_cross <- rbind(
  data.frame(go_inter_down_filtered, direction = "Hip < mPFC"),
  data.frame(go_inter_up_filtered, direction = "Hip > mPFC"))

#  data.frame(go_down_both_filtered, direction = "ctx_down_hp_down"),
#  data.frame(go_up_both_filtered, direction = "ctx_up_hp_up"))

go_combined_cross <- go_combined_cross %>%
  separate(
    GeneRatio, 
    into = c("GeneNum", "TotalNum"),  
    sep = "/", 
    convert = TRUE                   
  ) %>%
  mutate(
    GeneRatio = GeneNum / TotalNum    
  )

go_combined_cross_top15 <- go_combined_cross %>%
  group_by(direction) %>%
  arrange(p.adjust, .by_group = TRUE) %>%  
  slice_head(n = 20) %>%                    
  mutate(sort_priority = GeneRatio * Count) %>%
  arrange(desc(sort_priority), .by_group = TRUE) %>%  
  ungroup()

ggplot(go_combined_cross_top15, 
       aes(x = GeneRatio, 
           y = reorder(Description, sort_priority),  
           size = Count,
           color = -log10(p.adjust))) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(3, 8)) +
  facet_wrap(~ direction, scales = "free_y", ncol = 2) +  
  labs(
    title = "GO Enrichment of opposite KO effect in Cortex v.s. Hippocampus",
    x = "Gene Ratio", 
    y = "GO Term",
    size = "Gene Count",
    color = "-log10(padj)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )


# pheatmap

vsd <- vst(dds, blind = FALSE)
expr_matrix <- assay(vsd)
plot_matrix <- expr_matrix[rownames(expr_matrix) %in% union(rownames(shown_ctx), rownames(shown_hp)), ]

sample_info <- meta[,-1] # keep the "condition" & "region"
rownames(sample_info) <- colnames(plot_matrix) 

ann_colors = list(
  condition = c(WT = "#1B9E77", KO = "#D95F02"),
  region = c(ctx = "#7570B3", hp = "#E7298A")
)

colnames(plot_matrix) <- c("WT_mPFC1","WT_mPFC2","WT_mPFC3","WT_Hip1","WT_Hip2","WT_Hip3","KO_mPFC1","KO_mPFC2","KO_mPFC3","KO_Hip1","KO_Hip2","KO_Hip3")

target_order <- c("WT_mPFC1","WT_mPFC2","WT_mPFC3","KO_mPFC1","KO_mPFC2","KO_mPFC3","WT_Hip1","WT_Hip2","WT_Hip3","KO_Hip1","KO_Hip2","KO_Hip3")
plot_matrix_ordered <- plot_matrix[, target_order]
rownames(sample_info) <- colnames(plot_matrix) 

sample_info_ordered <- sample_info[target_order, ]

custom_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

pheatmap(plot_matrix_ordered, 
         scale = "row",             
         color = custom_colors,
         annotation_col = sample_info_ordered, 
         annotation_colors = ann_colors,
         cluster_cols = FALSE,      
         cluster_rows = TRUE,       
         show_rownames = TRUE,      
         show_colnames = TRUE,     
         fontsize_row = 10,
         border_color = "white",
         main = "Standardized Expression (Z-score) of Core Genes"
)

matrix_ctx <- plot_matrix_ordered[, 1:6]
matrix_hp <- plot_matrix_ordered[, 7:12]
sample_ctx <- sample_info_ordered[1:6, ]
sample_hp <- sample_info_ordered[7:12, ]

pheatmap(matrix_ctx, 
         scale = "row",             
         color = custom_colors,
         annotation_col = sample_ctx, 
         annotation_colors = ann_colors,
         cluster_cols = FALSE,      
         cluster_rows = TRUE,       
         show_rownames = TRUE,      
         show_colnames = TRUE,     
         fontsize_row = 10,
         border_color = "white",
         main = "Standardized Expression (Z-score) of Core Genes in mPFC"
)
pheatmap(matrix_hp, 
         scale = "row",             
         color = custom_colors,
         annotation_col = sample_hp, 
         annotation_colors = ann_colors,
         cluster_cols = FALSE,      
         cluster_rows = TRUE,       
         show_rownames = TRUE,      
         show_colnames = TRUE,     
         fontsize_row = 10,
         border_color = "white",
         main = "Standardized Expression (Z-score) of Core Genes in Hip"
)


# complexheatmap

matrix_mpfc_z <- t(scale(t(matrix_ctx)))
matrix_hip_z  <- t(scale(t(matrix_hp)))

mean_wt <- rowMeans(matrix_mpfc_z[, 1:3])
mean_ko <- rowMeans(matrix_mpfc_z[, 4:6])

gene_direction <- ifelse(mean_ko > mean_wt, "Shared Upregulated (KO > WT)", "Shared Downregulated (KO < WT)")
row_split_factor <- factor(gene_direction, levels = c("Shared Upregulated (KO > WT)", "Shared Downregulated (KO < WT)"))

anno_mpfc <- HeatmapAnnotation(
  Condition = factor(rep(c("WT", "KO"), each = 3), levels = c("WT", "KO")),
  Region = factor(rep("ctx", 6), levels = c("ctx", "hp")),
  col = list(
    Condition = c("WT" = "#1B9E77", "KO" = "#D95F02"), 
    Region = c("ctx" = "#7570B3", "hp" = "#E7298A")
  ),
  show_annotation_name = TRUE,          
  annotation_name_side = "left"
)

anno_hip <- HeatmapAnnotation(
  Condition = factor(rep(c("WT", "KO"), each = 3), levels = c("WT", "KO")),
  Region = factor(rep("hp", 6), levels = c("ctx", "hp")),
  col = list(
    Condition = c("WT" = "#1B9E77", "KO" = "#D95F02"), 
    Region = c("ctx" = "#7570B3", "hp" = "#E7298A")
  ),
  show_annotation_name = FALSE,         
  show_legend = FALSE                   
)

col_fun = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

ht_mpfc <- Heatmap(
  matrix_mpfc_z,
  name = "Z-score",
  col = col_fun,
  top_annotation = anno_mpfc,               
  
  cluster_columns = FALSE,
  cluster_rows = TRUE,                      
  row_split = row_split_factor,                 
  
  row_title = NULL,                         
  
  show_row_names = FALSE,                   
  
  show_column_names = TRUE,                 
  column_names_side = "top",               
  column_names_rot = 45,                    
  column_names_gp = gpar(fontsize = 10),
  
  rect_gp = gpar(col = "white", lwd = 1)
)

ht_hip <- Heatmap(
  matrix_hip_z,
  col = col_fun,
  top_annotation = anno_hip,                
  
  cluster_columns = FALSE,
  cluster_rows = TRUE,                      
  row_split = row_split_factor,                 
  row_title = NULL,                         
  
  show_heatmap_legend = FALSE,
  show_row_names = TRUE,                    
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 9, fontface = "italic"),
  show_column_names = TRUE,                 
  column_names_side = "top",
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10),
  
  rect_gp = gpar(col = "white", lwd = 1)
)

combined_ht_list <- ht_mpfc + ht_hip

draw(
  combined_ht_list, 
  main_heatmap = "Z-score",                 
  ht_gap = unit(2, "mm"),                   
  row_gap = unit(1, "mm")                 
)


# cnetplot 
both_genes <- res_df[res_df$Group %in% c("Both_up","Both_down"),]
saveRDS(both_genes,file = "both_genes.rds") 
both_genes$mean_log2FC <- rowMeans(both_genes[, c("lfc_ctx", "lfc_hp")])

fold_changes <- both_genes$mean_log2FC
names(fold_changes) <- both_genes$GeneSymbol

go_merged <- merge_result(
  list(
    Downregulated = go_shared_down,  
    Upregulated = go_shared_hp
  )
)

selected_pathways <- c("extracellular matrix organization", "cell-substrate adhesion", "regulation of interleukin-6 production",
                       "negative regulation of MAPK cascade","regulation of ERK1 and ERK2 cascade",
                       "regulation of axonogenesis", "regulation of neurogenesis", # Mbp
                       "fibroblast proliferation") # Esr1

gene_fc_status <- ifelse(fold_changes > 0, "Up", "Down")
gene_groups <- data.frame(
  gene = names(fold_changes),
  group = gene_fc_status,
  fc = fold_changes
)

p1 <- cnetplot(go_merged, 
               showCategory = selected_pathways, 
               circular = FALSE, 
               color.params = list(foldChange = fold_changes, edge = TRUE),
               repel = TRUE,              
               cex_label_category = 1.2,  
               cex_label_gene = 0.8       
)

p1 <- p1 + scale_color_gradient2(name = "Mean log2FC", 
                                 low = "blue",    
                                 mid = "white",   
                                 high = "red",    
                                 midpoint = 0)  

print(p1)
