#load libraries
library(stringr)
library(CellChat)
library(Seurat)
library(future)
library(dplyr)
library(readxl)
library(purrr)
library(tidyr)
library(openxlsx)
#######################load seurat data##################################
data_dir <- 'rawdata'
mat <- Read10X(data.dir = data_dir, gene.column = 1, cell.column = 1) 
meta_used <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
meta_used$Sample <- sapply(strsplit(rownames(meta_used), "\\."), `[`, 1)
snRNA <- CreateSeuratObject(counts = mat, project = "mdd", meta.data=meta_used)
colnames(snRNA) <- gsub('^"|"$', '', colnames(snRNA))
umap_data <- read.table(gzfile("UMAP.coords.tsv.gz"), row.names = 1)
colnames(umap_data) <- c("UMAP_1", "UMAP_2")
umap_reduc <- CreateDimReducObject(embeddings = as.matrix(umap_data), key = "UMAP_")
snRNA[["umap"]] <- umap_reduc
DimPlot(snRNA, reduction = "umap",group.by = "Cluster",label = TRUE)

saveRDS(snRNA,file = "snRNA.rds")
split_snRNA <- SplitObject(object = snRNA, split.by = "Sex")
names(split_snRNA) 
female_snRNA <- split_snRNA[["Female"]]
male_snRNA <- split_snRNA[["Male"]]
saveRDS(female_snRNA,file = "female_snRNA.rds")
saveRDS(male_snRNA,file = "male_snRNA.rds")


# ---- load Data ----
female_snRNA <- readRDS("female_snRNA.rds")
female_snRNA@meta.data$Sample <- str_extract(rownames(female_snRNA@meta.data), "^(F|M)\\d+")

male_snRNA <- readRDS("male_snRNA.rds")
male_snRNA@meta.data$Sample <- str_extract(rownames(male_snRNA@meta.data), "^(F|M)\\d+")

# ---- Set parallel processing ----
options(future.globals.maxSize = 27 * 1024^3)
#Linux or windows(multisession)
future::plan("multicore", workers = 8)

#############################Actual sample grouping####################################################
f_case_samps     <- c("F1", "F11", "F12", "F14", "F15", "F16", "F17", "F18", "F19",
                      "F2", "F20", "F25", "F27", "F28", "F3", "F4", "F5", "F6","F8","F9")
f_control_samps  <- c("F10", "F13", "F21", "F22", "F23", "F24", "F26", "F29",
                      "F30", "F31", "F32", "F33", "F34", "F35", "F36", "F37",
                      "F38", "F7")

m_case_samps <- c("M1", "M10", "M11", "M14", "M17", "M18", "M23", "M26", "M28", "M30", "M32", "M33", "M34", "M4", "M5", "M6", "M8")
m_control_samps <- c("M12", "M13", "M15", "M16", "M19", "M2", "M20", "M21", "M22", "M24", "M27", "M29", "M3", "M31", "M7", "M9")

# ---- Function: build_cellchat_from_samples ----
build_cellchat_from_samples <- function(seurat_obj, sample_ids) {
  cells_keep <- colnames(seurat_obj)[seurat_obj$Sample %in% sample_ids]
  subset_obj <- subset(seurat_obj, cells = cells_keep)
  Idents(subset_obj) <- "Cluster"
  
  cellchat <- createCellChat(object = subset_obj, group.by = "ident", assay = "RNA")
  CellChatDB <- CellChatDB.human
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use
  # CellChatDB.use <- subsetDB(cellchat@DB, search = list(c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), c("CellChatDB v1")), key = c("annotation", "version"))
  cellchat <- subsetData(cellchat)
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 5)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  return(cellchat)
}

#bulild CellChat object
f_ca_obj <- build_cellchat_from_samples(female_snRNA, f_case_samps)
f_ct_obj <- build_cellchat_from_samples(female_snRNA, f_control_samps)
m_ca_obj <- build_cellchat_from_samples(male_snRNA, m_case_samps)
m_ct_obj <- build_cellchat_from_samples(male_snRNA, m_control_samps)
  
#save the CellChat object as RDS
saveRDS(f_ca_obj,file = "Female_MDD.rds")
saveRDS(f_ct_obj,file = "Female_control.rds") 
saveRDS(m_ca_obj,file = "Male_MDD.rds")
saveRDS(m_ct_obj,file = "Male_control.rds") 



####################cellchat_list
snRNA <- readRDS("snRNA.rds")
snRNA@meta.data$Sample <- str_extract(rownames(snRNA@meta.data), "^(F|M)\\d+")
Idents(object = snRNA) <- "Cluster"
table(snRNA@meta.data$Sample)
sample_list <- list()
unique_samples <- unique(snRNA@meta.data$Sample)
for (sample in unique_samples) {
  sample_cells <- rownames(snRNA@meta.data[snRNA@meta.data$Sample == sample, ])
  sample_obj <- subset(snRNA, cells = sample_cells)
  sample_list[[sample]] <- sample_obj
}
cellchat_list <- list()
for (sample in names(sample_list)) {
  seurat_obj <- sample_list[[sample]]
  set.seed(12345)
  data.input <- seurat_obj[["RNA"]]@data
  labels <- Idents(seurat_obj)
  meta <- data.frame(labels = labels, row.names = names(labels))
  cellchat <- createCellChat(object = seurat_obj, group.by = "Cluster", assay = "RNA")
  CellChatDB <- CellChatDB.human
  showDatabaseCategory(CellChatDB)
  CellChatDB.use <- subsetDB(CellChatDB)
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  future::plan("multisession", workers = 12)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  ptm <- Sys.time()
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  execution.time <- Sys.time() - ptm
  print(paste("Execution time for", sample, ":", as.numeric(execution.time, units = "secs"), "seconds"))
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  cellchat_list[[sample]] <- cellchat
}

saveRDS(cellchat_list,file = "cellchat_list.rds")










