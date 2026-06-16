##################################permutationTest##########################################################


args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("USECASE: Rscript permutationTest_modified.R <start_id> <end_id> <combinations_file>")
}

start_id <- as.integer(args[1])
end_id <- as.integer(args[2])
combinations_file <- args[3]


library(dplyr)
library(NeuronChat)
library(CellChat)
library(Seurat)
library(ggplot2)
library(ggalluvial)
library(ComplexHeatmap)
library(circlize)
library(scales)

# ---- Set working directory and load data ----
# setwd("/home/liyinong/permutationTest/neuronChatF")
setwd("neuronChatF")
source('heatmap_aggregated_v2.R')
source('net_aggregation.R')

cat("Reading SampleData:", combinations_file, "\n")
combinations_data <- read.csv(combinations_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8")

cat("success Reading SampleData", nrow(combinations_data), " datas\n")
cat("the range of current task: ID", start_id, "to", end_id, "\n")
# group info
cluster_to_group <- data.frame(
  Cluster = c("Ast1", "Ast2", "End1", "ExN1_L24","ExN10_L46","ExN11_L56","ExN12_L56","ExN13_L56","ExN14","ExN15_L56","ExN16_L56","ExN17","ExN18","ExN19_L56","ExN2_L23","ExN20_L56","ExN3_L46","ExN4_L35","ExN5","ExN6","ExN7", "ExN8_L24","ExN9_L23","InN1_PV","InN10_ADARB2","InN2_SST","InN3_VIP","InN4_VIP","InN5_SST","InN6_LAMP5","InN7_Mix","InN8_Mix","InN9_PV","Mic1","Mix","Oli1","Oli2","Oli3","OPC1","OPC2","OPC3"),
  Group = c("Non-Neuronal", "Non-Neuronal", "Non-Neuronal", "Glutamatergic", "Glutamatergic", "Glutamatergic", "Glutamatergic","Glutamatergic", "Glutamatergic", "Glutamatergic", "Glutamatergic","Glutamatergic", "Glutamatergic", "Glutamatergic", "Glutamatergic","Glutamatergic", "Glutamatergic", "Glutamatergic", "Glutamatergic","Glutamatergic", "Glutamatergic", "Glutamatergic", "Glutamatergic","GABAergic","GABAergic","GABAergic","GABAergic","GABAergic","GABAergic","GABAergic","GABAergic","GABAergic","GABAergic","Non-Neuronal", "Non-Neuronal", "Non-Neuronal","Non-Neuronal", "Non-Neuronal", "Non-Neuronal","Non-Neuronal", "Non-Neuronal")
)
group <- setNames(cluster_to_group$Group, cluster_to_group$Cluster)

############################## OR reading original data ##################################
data_dir <- 'rawdata'
mat <- Read10X(data.dir = data_dir, gene.column = 1, cell.column = 1) 
meta_used <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
meta_used$Sample <- sapply(strsplit(rownames(meta_used), "\\."), `[`, 1)

meta <- meta_used %>% left_join(cluster_to_group, by= "Cluster")
rownames(meta) <- rownames(meta_used)

# ---- get case and control from SampleData ----
get_samples_from_combination <- function(combinations_data, comb_id) {
  if (comb_id > nrow(combinations_data)) {
    stop("The requested ID is invalid and out of range :", comb_id)
  }
  
  case_samples <- trimws(unlist(strsplit(combinations_data$case_samples[comb_id], ",")))
  control_samples <- trimws(unlist(strsplit(combinations_data$control_samples[comb_id], ",")))
  
  return(list(case = case_samples, control = control_samples))
}

######################################Permutation cycle ##################################################
# ---- initialization parameter ----

cat("Start performing the permutation test.....\n")

for (comb_id in start_id:end_id) {
  if (comb_id > nrow(combinations_data)) {
    cat("Warning: The combination ID", comb_id, "is out of the data range.Skiping\n")
    next
  }
  
  cat("Dealing with the combination", comb_id, "/", nrow(combinations_data), "\n")
  iter_dir_name <- paste0("iter_", comb_id, "_", gsub("[: -]", "", Sys.time()), "_", Sys.getpid())
  dir.create(iter_dir_name, recursive = TRUE, showWarnings = FALSE)
  setwd(iter_dir_name)
  
  if ((comb_id - start_id + 1) %% 10 == 0) {
    cat("Progress:", (comb_id - start_id + 1), "/", (end_id - start_id + 1), 
        "Total:", comb_id, "/", nrow(combinations_data), "\n")
  }
  
  # 1. Permutation condition label
  samples <- get_samples_from_combination(combinations_data, comb_id)
  f_case_samps <- samples$case
  f_control_samps <- samples$control
  record_data <- data.frame(
    combination_id = comb_id,
    case_samples = paste(f_case_samps, collapse = ","),
    control_samples = paste(f_control_samps, collapse = ","),
    timestamp = Sys.time()
  )
  write.csv(record_data, file = "combination_detail_record.csv", row.names = FALSE)
  
  # 2. bulild NeuronChat object
  cell_fca_in <- which(meta$Sample %in% f_case_samps)
  cell_fct_in <- which(meta$Sample %in% f_control_samps)
  meta_fca <- meta[cell_fca_in,]
  meta_fct <- meta[cell_fct_in,]
  
  target_fca <- mat[, cell_fca_in]
  target_fct <- mat[, cell_fct_in]
  
  set.seed(1234)
  fca_x <- createNeuronChat(as.matrix(target_fca), DB='human',group.by = meta_fca$Cluster,meta=meta_fca);
  fca_x <- run_NeuronChat(fca_x,M=10)
  set.seed(1234)
  fct_x <- createNeuronChat(as.matrix(target_fct), DB='human',group.by = meta_fct$Cluster,meta=meta_fct);
  fct_x <- run_NeuronChat(fct_x,M=10)
  
  rm(meta_fca,meta_fct,target_fct,target_fca)
  gc()
  
  # 3. data processing
  data_list <- c("Female_Control","Female_MDD")
  nc <- mergeNeuronChat(list(fct_x,fca_x), add.names = data_list)
  
  # global heatmap for female
  heatmap_aggregated_v2(nc, dataset=data_list,method='weight',group = group)
  
  #################### with NRXN3_NLGN1 or NRXN1_NLGN1 only ################

  fca_n1 <- fca_x@net[["NRXN1_NLGN1"]]
  fct_n1 <- fct_x@net[["NRXN1_NLGN1"]]
  fca_n3 <- fca_x@net[["NRXN3_NLGN1"]]
  fct_n3 <- fct_x@net[["NRXN3_NLGN1"]]
  
  nn <- list(fct_n3,fca_n3) # for NRXN3_NLGN1
  nn <- list(fct_n1,fca_n1) # for NRXN1_NLGN1

  val <- data.frame(Data=numeric(), Cluster= numeric(), Value=numeric(),stringsAsFactors = F)
  for (i in 1: length(nn)){   # calculating female datasets
    print(i)
    sender <- nn[[i]][36:41,]
    receiver <- nn[[i]][,36:41]
    
    print(rowSums(sender))
    print(colSums(receiver))
    tmp <- data.frame(Data=i,Cluster=names(rowSums(sender)), Value=rowSums(sender))
    val <- rbind(val,tmp)
    tmp <- data.frame(Data=i,Cluster=names(colSums(receiver)), Value=colSums(receiver))
    val <- rbind(val,tmp)
  }
  
  
  n_datasets <- 2
  n_sexes <- 1
  n_celltypes <- 6
  
  
  data <- expand.grid(
    CellType = c("OPC1","OPC2","OPC3","Oli1","Oli2","Oli3"),
    Direction = c("Outgoing", "Incoming"),
    Dataset = c("Control","MDD")
  ) %>%
    mutate(
      CellType = factor(CellType, levels = unique(CellType)),
      Dataset = factor(Dataset, levels = c("Control", "MDD")),
      Value = val$Value)
  
  data <- data %>% mutate(
    Value = ifelse(Direction == "Outgoing", -1 * abs(Value), abs(Value)))
  
  data <- data %>% mutate(
    GroupID = interaction(
      CellType, Direction, # CellType
      Dataset,
      sep = "_",
      lex.order = FALSE 
    )
  ) 
  
  
  data <- data %>% mutate(
    Show = interaction(Dataset,sep = "_")
  )
  
  head(data)
  
  
  celltype_spacing <- 1  
  dataset_offset <- 0.2   
  
  data <- data %>%
    mutate(
      x_base = as.numeric(CellType) * celltype_spacing,
    )
  
  data <- data %>%
    group_by(CellType) %>%
    mutate(
      x_offset = (as.numeric(Show) - (n_datasets+1)/2) * dataset_offset,
      x_final = x_base + x_offset
    ) 
  
  
  write.csv(data, file = "data.csv", row.names = FALSE)
  rm(fca_x,fct_x,fca_n1,fct_n1,fca_n3,fct_n3,nn,nc)
  setwd("..")
  gc() 
  
  
}
cat("Task completed! Processed combinations", start_id, "to", end_id, "\n")












