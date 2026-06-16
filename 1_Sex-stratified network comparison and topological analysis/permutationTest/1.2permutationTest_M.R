##################################permutationTest##########################################################


args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("USECASE: Rscript permutationTest_modified.R <start_id> <end_id> <combinations_file>")
}

start_id <- as.integer(args[1])
end_id <- as.integer(args[2])
combinations_file <- args[3]

library(stringr)
library(CellChat)
library(Seurat)
library(future)
library(dplyr)
library(readxl)
library(purrr)
library(tidyr)
library(openxlsx)

# ---- Set working directory and load data ----
setwd("/home/liyinong/permutationTest/male")
source('ChangesScatter-for permutation.R')
source('contribution-for permutation.R')
source('scatter-for permutation.R')
# ---- load SampleData ----
cat("Reading SampleData:", combinations_file, "\n")
combinations_data <- read.csv(combinations_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8")

cat("success Reading SampleData", nrow(combinations_data), " datas\n")
cat("the range of current task: ID", start_id, "to", end_id, "\n")


male_snRNA <- readRDS("male_snRNA.rds")
male_snRNA@meta.data$Sample <- str_extract(rownames(male_snRNA@meta.data), "^(F|M)\\d+")

# ---- Set parallel processing ----
options(future.globals.maxSize = 27 * 1024^3)
future::plan("multicore", workers = 8)

#############################Actual sample grouping####################################################
# male_case_list <- c("M1", "M10", "M11", "M14", "M17", "M18", "M23", "M26", "M28", "M30", "M32", "M33", "M34", "M4", "M5", "M6", "M8")
# male_control_list <- c("M12", "M13", "M15", "M16", "M19", "M2", "M20", "M21", "M22", "M24", "M27", "M29", "M3", "M31", "M7", "M9")

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

# ---- Function: attenuation analysis  ----

compute_global_strength <- function(probability_matrix) { sum(abs(probability_matrix)) }

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
  m_case_samps <- samples$case
  m_control_samps <- samples$control
  record_data <- data.frame(
    combination_id = comb_id,
    case_samples = paste(m_case_samps, collapse = ","),
    control_samples = paste(m_control_samps, collapse = ","),
    timestamp = Sys.time()
  )
  write.csv(record_data, file = "combination_detail_record.csv", row.names = FALSE)
  # 2. bulild CellChat object
  m_ca_obj <- build_cellchat_from_samples(male_snRNA, m_case_samps)
  m_ct_obj <- build_cellchat_from_samples(male_snRNA, m_control_samps)
  
  saveRDS(m_ca_obj,
        file = "male_case_cellchat.rds")


  saveRDS(m_ct_obj,
        file = "male_control_cellchat.rds")  
  
  # 3. attenuation analysis
  m_ca_mat <- m_ca_obj@net$prob
  m_ct_mat <- m_ct_obj@net$prob
  
  global_male_control <- compute_global_strength(m_ct_mat)
  global_male_mdd <- compute_global_strength(m_ca_mat)
  
  male_attenuation <- ifelse(global_male_control != 0, (1 - global_male_mdd / global_male_control) * 100, 0)
  #########################result1
  message("####################male_attenuation####################",male_attenuation)
  attenuation <- data.frame(
    combination_id = comb_id,
    attenuationValue = male_attenuation
  )
  write.csv(attenuation, file = "attenuation.csv", row.names = FALSE)
  
  
  ############################################step2####################################################  
  object.list <- list(male_control = m_ct_obj, male_MDD = m_ca_obj)
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))
  # rm(object.list)
  
  # ref: scatter-for permutation.R
  netAnalysis_signalingRole_scatter_V2(cellchat, data=c("male_control","male_MDD"), F_flag = F, M_flag= T,external_objects = object.list)

  ############################################step5####################################################  
  male_case_pathways <- m_ca_obj@netP$pathways
  male_control_pathways <- m_ct_obj@netP$pathways
  all_pathways <- unique(c(male_case_pathways, male_control_pathways))
  calculate_pathway_sums <- function(obj) {
    pathways <- obj@netP$pathways
    prob_array <- obj@netP$prob
    sums <- numeric(length(pathways))
    names(sums) <- pathways
    
    for (i in seq_along(pathways)) {
      pathway <- pathways[i]
      mat <- prob_array[, , pathway]
      sums[i] <- sum(mat)
    }
    
    return(sums)
  }
  male_case_sums <- calculate_pathway_sums(m_ca_obj)
  male_control_sums <- calculate_pathway_sums(m_ct_obj)
  
  result_df <- data.frame(pathway = all_pathways)
  
  
  result_df$male_case <- sapply(all_pathways, function(p) {
    if (p %in% names(male_case_sums)) male_case_sums[p] else 0
  })
  
  result_df$male_control <- sapply(all_pathways, function(p) {
    if (p %in% names(male_control_sums)) male_control_sums[p] else 0
  })
  
  
  result_df$male_diff <- 1 - (result_df$male_case / result_df$male_control)
  result_df$male_diff[!is.finite(result_df$male_diff)] <- NA
  result_df$caseMinuscontrol <- result_df$male_case - result_df$male_control
  print(head(result_df))
  write.csv(result_df, "venn_dif_3.csv", row.names = FALSE)  
  
  rm(m_ct_obj,m_ca_obj)
  gc()
  
  file <- "Differences.xlsx"
  sheet_names <- excel_sheets(file)
  print(sheet_names)  
  
  data_list <- map(
    sheet_names, ~ read_excel( path = file, sheet = .x, col_types = "text")
  ) 
  
  names(data_list) <- sheet_names
  data_list <- data_list[5:7]  # sheet1: overall; sheet2-4: ed/x/y; sheet5-7: ed/x/y top10
  
  #male_cell_types <- c("Ast1", "Ast2", "End1", "ExN1_L24", "ExN10_L46", "ExN11_L56", "ExN12_L56", "ExN13_L56", "ExN14", "ExN15_L56", "ExN16_L56", "ExN18", "ExN19_L56", "ExN2_L23", "ExN20_L56", "ExN3_L46", "ExN4_L35", "ExN5", "ExN6", "ExN7", "ExN8_L24", "ExN9_L23", "InN1_PV", "InN10_ADARB2", "InN2_SST", "InN3_VIP", "InN4_VIP", "InN5_SST", "InN6_LAMP5", "InN7_Mix", "InN8_Mix", "InN9_PV", "Mic1", "Mix", "Oli1", "Oli2", "Oli3", "OPC1", "OPC2", "OPC3")  
  #female_cell_types <- c("Ast1", "Ast2", "End1", "ExN1_L24", "ExN10_L46", "ExN11_L56", "ExN12_L56", "ExN13_L56", "ExN14", "ExN15_L56", "ExN16_L56", "ExN17", "ExN18", "ExN19_L56", "ExN2_L23", "ExN20_L56", "ExN3_L46", "ExN4_L35", "ExN5", "ExN6", "ExN7", "ExN8_L24", "ExN9_L23", "InN1_PV", "InN10_ADARB2", "InN2_SST", "InN3_VIP", "InN4_VIP", "InN5_SST", "InN6_LAMP5", "InN7_Mix", "InN8_Mix", "InN9_PV", "Mic1", "Mix", "Oli1", "Oli2", "Oli3", "OPC1", "OPC2", "OPC3")  
  
  process_sheet <- function(sheet_data, direction_type) {
    sheet_data %>%
      mutate(cell_type = factor(cell_type, levels = cell_type),
             sth = as.numeric(sth),
             direction = direction_type) %>%
      select(cell_type, direction)
  }
  
  processed_list <- list()
  
  for(sheet in names(data_list)){
    meta <- str_split(sheet, "_", simplify = TRUE)
    current_direction <- meta[1]
    
    processed <- process_sheet(
      data_list[[sheet]],
      direction_type = current_direction
    )
    print(processed)
    processed_list[[sheet]] <- processed
  }
  
  combined_data <- bind_rows(processed_list)
  
  combined_data <- combined_data %>%
    mutate(
      cell_mapping = case_when(
        str_detect(cell_type, "Ast") ~ "Ast" ,
        str_detect(cell_type, "End") ~ "End" ,
        str_detect(cell_type, "ExN") ~ "ExN" ,
        str_detect(cell_type, "InN") ~ "InN" ,
        str_detect(cell_type, "Mic") ~ "Mic" ,
        str_detect(cell_type, "Mix") ~ "Mix" ,
        str_detect(cell_type, "OPC") ~ "OPC" ,
        str_detect(cell_type, "Oli") ~ "OL" ,
        TRUE ~ "Other" 
      ))  
  
  
  ##################################### F2a2b2c Ratio
  
  count_opc_ol <- function(df) {
    opc_count <- sum(grepl("OPC", df$cell_mapping, ignore.case = TRUE))
    ol_count <- sum(grepl("OL", df$cell_mapping, ignore.case = TRUE))
    return(c(OPC = opc_count, OL = ol_count))
  }
  
  layer_counts <- combined_data %>%
    group_by(direction) %>%
    group_modify(~ as.data.frame(t(count_opc_ol(.x))))
  
  wide_counts <- layer_counts %>%
    pivot_wider(
      names_from = direction,
      values_from = c(OPC, OL)
    )
  
  Counts <- c(
    wide_counts$OPC_Euclidean,  
    wide_counts$OL_Euclidean,   
    wide_counts$OPC_Horizontal,   
    wide_counts$OL_Horizontal,    
    wide_counts$OPC_Vertical,   
    wide_counts$OL_Vertical     
  )
  
  names(Counts) <- c("OPC_ed", "OL_ed", "OPC_x", "OL_x", "OPC_y", "OL_y")
  Counts_df <- data.frame(
    Category = names(Counts),
    Value = as.numeric(Counts)
  )
  write.csv(Counts_df, file = "Counts_results.csv", row.names = FALSE)
  # original counts: c(1,1,1,0,2,2)
  Total <- rep(c(3,3), 3)                          # No. of superclass (fixed)
  Ratio <- round(Counts / Total,2)
  
  Counts <- rbind(Counts, Total, Ratio)
  
  print(Counts)
  
  ############################################step3####################################################   
  sheet_names <- excel_sheets(file)
  print(sheet_names)  
  
  data_list <- map(
    sheet_names, ~ read_excel( path = file, sheet = .x, col_types = "text")
  ) 
  
  names(data_list) <- sheet_names
  data_list <- data_list[5]  # sheet1: overall; sheet2-4: ed/x/y; sheet5-7: ed/x/y top10
  
  # parsing sheet5 for ed_top_10 clusters
  sheet <- names(data_list[[1]])
  sheet_data <- data_list[[1]]
  #sheet_data <- sheet_data %>% mutate(cell_type = factor(cell_type, levels = cell_type)) 
  print(sheet_data$cell_type)
  netAnalysis_signalingChanges_scatter_V2(cellchat,idents.use.multi = sheet_data$cell_type)
  #c("ExN7", "OPC1", "ExN17", "InN3_VIP", "ExN8_L24", "ExN14", "ExN13_L56", "Oli1", "ExN18", "ExN11_L56"))
  
  ###################################### F2h2i signalingChanges_scatter
  
  netAnalysis_contribution_v2(cellchat, dataset = c("male_control","male_MDD"), signaling = "NRXN")
  setwd("..")
  gc()
}
cat("Task completed! Processed combinations", start_id, "to", end_id, "\n")