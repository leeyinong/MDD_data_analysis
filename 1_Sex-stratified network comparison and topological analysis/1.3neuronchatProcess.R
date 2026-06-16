#load libraries
library(dplyr)
library(NeuronChat)
library(CellChat)
library(Seurat)
library(ggplot2)
library(ggalluvial)
library(ComplexHeatmap)
library(circlize)

#load Data
data_dir <- 'rawdata'
mat <- Read10X(data.dir = data_dir, gene.column = 1, cell.column = 1) 
meta_used <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
meta_used$Sample <- sapply(strsplit(rownames(meta_used), "\\."), `[`, 1)

# group info
cluster_to_group <- data.frame(
  Cluster = c("Ast1", "Ast2", "End1", "ExN1_L24","ExN10_L46","ExN11_L56","ExN12_L56","ExN13_L56","ExN14","ExN15_L56","ExN16_L56","ExN17","ExN18","ExN19_L56","ExN2_L23","ExN20_L56","ExN3_L46","ExN4_L35","ExN5","ExN6","ExN7", "ExN8_L24","ExN9_L23","InN1_PV","InN10_ADARB2","InN2_SST","InN3_VIP","InN4_VIP","InN5_SST","InN6_LAMP5","InN7_Mix","InN8_Mix","InN9_PV","Mic1","Mix","Oli1","Oli2","Oli3","OPC1","OPC2","OPC3"),
  Group = c("Non-Neuronal", "Non-Neuronal", "Non-Neuronal", "Glutamatergic", "Glutamatergic", "Glutamatergic", "Glutamatergic","Glutamatergic", "Glutamatergic", "Glutamatergic", "Glutamatergic","Glutamatergic", "Glutamatergic", "Glutamatergic", "Glutamatergic","Glutamatergic", "Glutamatergic", "Glutamatergic", "Glutamatergic","Glutamatergic", "Glutamatergic", "Glutamatergic", "Glutamatergic","GABAergic","GABAergic","GABAergic","GABAergic","GABAergic","GABAergic","GABAergic","GABAergic","GABAergic","GABAergic","Non-Neuronal", "Non-Neuronal", "Non-Neuronal","Non-Neuronal", "Non-Neuronal", "Non-Neuronal","Non-Neuronal", "Non-Neuronal")
)
group <- setNames(cluster_to_group$Group, cluster_to_group$Cluster)

meta <- meta_used %>% left_join(cluster_to_group, by= "Cluster")
rownames(meta) <- rownames(meta_used)

# bulild NeuronChat object
cell_fca_in <- which(meta$Sex=="Female" & meta$Condition=="Case")
cell_fct_in <- which(meta$Sex=="Female" & meta$Condition=="Control")

cell_mca_in <- which(meta$Sex=="Male" & meta$Condition=="Case")
cell_mct_in <- which(meta$Sex=="Male" & meta$Condition=="Control")

meta_fca <- meta[cell_fca_in,]
meta_fct <- meta[cell_fct_in,]

meta_mca <- meta[cell_mca_in,]
meta_mct <- meta[cell_mct_in,]

target_fca <- mat[, cell_fca_in]
target_fct <- mat[, cell_fct_in]

target_mca <- mat[, cell_mca_in]
target_mct <- mat[, cell_mct_in]

set.seed(1234)
mca_x <- createNeuronChat(as.matrix(target_mca), DB='human',group.by = meta_mca$Cluster,meta=meta_mca);
mca_x <- run_NeuronChat(mca_x,M=10)
saveRDS(mca_x,file = "male_case_neuronchat.rds")
set.seed(1234)
mct_x <- createNeuronChat(as.matrix(target_mct), DB='human',group.by = meta_mct$Cluster,meta=meta_mct);
mct_x <- run_NeuronChat(mct_x,M=10)
saveRDS(mct_x,file = "male_control_neuronchat.rds")

set.seed(1234)
fca_x <- createNeuronChat(as.matrix(target_fca), DB='human',group.by = meta_fca$Cluster,meta=meta_fca);
fca_x <- run_NeuronChat(fca_x,M=10)
saveRDS(fca_x,file = "female_case_neuronchat.rds")
set.seed(1234)
fct_x <- createNeuronChat(as.matrix(target_fct), DB='human',group.by = meta_fct$Cluster,meta=meta_fct);
fct_x <- run_NeuronChat(fct_x,M=10)
saveRDS(fct_x,file = "female_control_neuronchat.rds")

