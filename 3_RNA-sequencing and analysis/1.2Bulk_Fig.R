#load libraries
library(tidyverse)
library(ggplot2)
library(patchwork)


#load Data
up   <- readxl::read_xlsx("go_shared_hp.xlsx")
down <- readxl::read_xlsx("go_shared_down.xlsx")

down_top30 <- down %>% arrange(desc(Count)) %>% slice_head(n = 30)
up_top30   <- up %>% arrange(desc(Count)) %>% slice_head(n = 30)

down_top30 <- down_top30 %>%
  arrange(Count, desc(p.adjust)) %>%
  mutate(Description = factor(Description, levels = Description))

up_top30 <- up_top30 %>%
  arrange(Count, desc(p.adjust)) %>%
  mutate(Description = factor(Description, levels = Description))


max_count <- max(c(up_top30$Count, down_top30$Count))
x_breaks  <- seq(0, ceiling(max_count), by = 2)
x_breaks_left  <- seq(0, 7, by = 2)
p_left <- ggplot(down_top30, aes(x = -Count, y = Description)) +
  geom_col(aes(fill = -log10(p.adjust)), width = 0.7, color = "white") +
  geom_text(
    aes(label = Description),x = -0.3,hjust = 1,size = 3,color = "black"
  ) +
  scale_fill_gradient(low = "#E6F0FA", high = "#4682B4", name = "-log10(P)", labels = function(x) sprintf("%.1f", x)) +
  labs(x = "Gene Count", y = NULL) +
  scale_x_continuous(limits = c(-7, 0), breaks = -x_breaks_left, labels = abs, expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.line.x = element_line(linewidth = 0.5), 
    axis.ticks.x = element_line(linewidth = 0.5), 
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(r = 0, l = 0, t = 5, b = 5)
  ) +
  ggtitle("Down-regulated")

p_right <- ggplot(up_top30, aes(x = Count, y = Description)) +
  geom_col(aes(fill = -log10(p.adjust)), width = 0.7, color = "white") +
  geom_text(
    aes(label = Description),x = 0.3,hjust = 0, size = 3,color = "black"
  ) +
  scale_fill_gradient(low = "#FFE6E6", high = "#D62728", name = "-log10(P)", labels = function(x) sprintf("%.1f", x)) +
  labs(x = "Gene Count", y = NULL) +
  scale_x_continuous(labels = abs,expand = c(0, 0),limits = c(0, max_count), breaks = x_breaks,) + 
  theme_minimal() +
  theme(
    axis.line.x = element_line(linewidth = 0.5),
    axis.ticks.x = element_line(linewidth = 0.5), 
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  ggtitle("Up-regulated")
  
p_combined <- p_left + p_right +
  plot_layout(widths = c(1, 1.58),guides = "collect") &
  theme(legend.position = "bottom",panel.spacing.x = unit(0, "pt"))

p_combined
#####################################################
library(ggtangle)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
library(scales)
# cnetplot 
both_genes <- readRDS("both_genes.rds")
go_shared_down <- readRDS("go_shared_down.rds")
go_shared_hp <- readRDS("go_shared_hp.rds")

fold_changes <- both_genes$mean_log2FC
names(fold_changes) <- both_genes$gene

go_merged <- merge_result(
  list(
    Upregulated = go_shared_hp,
    Downregulated = go_shared_down  
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

selected_pathways <- intersect(selected_pathways, go_merged@compareClusterResult$Description)


p1 <- cnetplot(go_merged, 
               showCategory = selected_pathways, 
               foldChange = fold_changes
)

df_pvalue <- go_merged@compareClusterResult[, c("Cluster", "Description", "p.adjust")] %>%
  filter(Description %in% selected_pathways)

# trans to *
df_pvalue$stars <- sapply(df_pvalue$p.adjust, function(p) {
  if (p < 0.0001) return("****")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
})

node_coords <- p1$data[1:8, c("name", "x", "y", "size")]
df_stars_plot <- merge(df_pvalue, node_coords, by.x = "Description", by.y = "name")

df_stars_plot <- df_stars_plot %>%
  mutate(
    x_offset = ifelse(Cluster == "Upregulated", x - 0.15, x + 0.15),
    y_offset = y + 0.05 
  )
print(p1)

display_min <- -2
display_max <- 3
p_fixed <- p1 + 
  scale_color_gradientn(
    name = "Fold Change",
    colors = c("#4682B4", "white", "#D62728"),
    limits = c(display_min, display_max),
    values = rescale(c(display_min, 0, display_max)),
    na.value = "#B3B3B3",
    oob = scales::squish,
  ) +
  guides(alpha = "none", 
         color = guide_colorbar(order = 1)) +
  geom_text(data = df_stars_plot, 
            aes(x = x_offset, y = y_offset, label = stars, group = Cluster), 
            color = "black", 
            fontface = "bold", 
            size = 4)

print(p_fixed)


