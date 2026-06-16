library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(openxlsx)


setwd("female")
# 1. Get all result folders
iter_folders <- list.dirs(path = ".", recursive = FALSE, full.names = FALSE)
iter_folders <- iter_folders[grepl("^iter_\\d+_", iter_folders)]
cat("find", length(iter_folders), "result folders\n")

#Check whether the result file is complete
check_files <- function(folders, required_files = NULL) {
  if (is.null(required_files)) {
    required_files <- c("attenuation.csv", "Counts_results.csv", 
                        "Pathways.xlsx", "Pairs.xlsx")
  }
  missing_info <- list()
  
  for (folder in folders) {
    missing_in_folder <- c()
    for (file in required_files) {
      file_path <- file.path(folder, file)
      if (!file.exists(file_path)) {
        missing_in_folder <- c(missing_in_folder, file)
      }
    }
    if (length(missing_in_folder) > 0) {
      missing_info[[folder]] <- missing_in_folder
    }
  }
  return(missing_info)
}
missing_files_list <- check_files(iter_folders)


# Initialize the summary table
attenuation_summary <- data.frame()
counts_summary <- data.frame()
nrxn_summary <- data.frame()
others_summary <- data.frame()
nrxn3_nlgn1_summary <- data.frame()
nrxn1_nlgn1_summary <- data.frame()
others_pairs_summary <- data.frame()


# 2. Traverse each folder
for (i in seq_along(iter_folders)) {
  folder <- iter_folders[i]
  cat("Dealing with folder:", folder, "...\n")
  
  # extract comb_id
  comb_id <- gsub("^iter_(\\d+)_.*", "\\1", folder)
  
  # file path
  attenuation_file <- file.path(folder, "attenuation.csv")
  counts_file <- file.path(folder, "Counts_results.csv")
  pathways_file <- file.path(folder, "Pathways.xlsx")
  pairs_file <- file.path(folder, "Pairs.xlsx")
  
  # 2.1 deal with attenuation.csv
  if (file.exists(attenuation_file)) {
    attenuation_data <- read.csv(attenuation_file)
    if (i == 1) {
      attenuation_summary <- attenuation_data
    } else {
      attenuation_summary <- bind_rows(attenuation_summary, attenuation_data)
    }
  }
  
  # 2.2 deal with Counts_results.csv
  if (file.exists(counts_file)) {
    counts_data <- read.csv(counts_file)
    
    # Transpose the data and convert the "Value" column into a column with the name "comb_id"
    counts_wide <- counts_data %>%
      select(Category, Value) %>%
      spread(key = Category, value = Value)
    
    # Add the "comb_id" column
    counts_wide <- counts_wide %>%
      mutate(comb_id = comb_id) %>%
      select(comb_id, everything())
    
    if (i == 1) {
      counts_summary <- counts_wide
    } else {
      counts_summary <- bind_rows(counts_summary, counts_wide)
    }
  }
  
  # 2.3 deal with Pathways.xlsx
  if (file.exists(pathways_file)) {
    # deal with NRXN sheet
    nrxn_data <- tryCatch({
      read_excel(pathways_file, sheet = "NRXN")
    }, error = function(e) {
      cat("folder", folder, "Not found NRXN sheet\n")
      return(NULL)
    })
    
    if (!is.null(nrxn_data) && nrow(nrxn_data) > 0) {
      # Keep only the required columns
      nrxn_simple <- nrxn_data %>%
        select(pathway, ed) %>%
        mutate(comb_id = comb_id) %>%
        select(comb_id, pathway, ed)
      
      if (i == 1) {
        nrxn_summary <- nrxn_simple
      } else {
        nrxn_summary <- bind_rows(nrxn_summary, nrxn_simple)
      }
    }
    
    # deal with Others sheet
    others_data <- tryCatch({
      read_excel(pathways_file, sheet = "Others")
    }, error = function(e) {
      cat("folder", folder, "Not found Others sheet\n")
      return(NULL)
    })
    
    if (!is.null(others_data) && nrow(others_data) > 0) {
      # Keep only the required columns
      others_simple <- others_data %>%
        select(pathway, ed) %>%
        mutate(comb_id = comb_id) %>%
        select(comb_id, pathway, ed)
      
      if (i == 1) {
        others_summary <- others_simple
      } else {
        others_summary <- bind_rows(others_summary, others_simple)
      }
    }
  }
  
  # 2.4 deal with Pairs.xlsx
  if (file.exists(pairs_file)) {
    # deal with NRXN3_NLGN1 sheet
    nrxn3_nlgn1_data <- tryCatch({
      read_excel(pairs_file, sheet = "NRXN3_NLGN1")
    }, error = function(e) {
      cat("folder", folder, "Not found NRXN3_NLGN1 sheet\n")
      return(NULL)
    })
    
    if (!is.null(nrxn3_nlgn1_data) && nrow(nrxn3_nlgn1_data) > 0) {
      nrxn3_nlgn1_simple <- nrxn3_nlgn1_data %>%
        select(name, contribution) %>%
        mutate(comb_id = comb_id) %>%
        select(comb_id, name, contribution)
      
      if (i == 1) {
        nrxn3_nlgn1_summary <- nrxn3_nlgn1_simple
      } else {
        nrxn3_nlgn1_summary <- bind_rows(nrxn3_nlgn1_summary, nrxn3_nlgn1_simple)
      }
    }
    
    # deal with NRXN1_NLGN1 sheet
    nrxn1_nlgn1_data <- tryCatch({
      read_excel(pairs_file, sheet = "NRXN1_NLGN1")
    }, error = function(e) {
      cat("folder", folder, "Not found NRXN1_NLGN1 sheet\n")
      return(NULL)
    })
    
    if (!is.null(nrxn1_nlgn1_data) && nrow(nrxn1_nlgn1_data) > 0) {
      nrxn1_nlgn1_simple <- nrxn1_nlgn1_data %>%
        select(name, contribution) %>%
        mutate(comb_id = comb_id) %>%
        select(comb_id, name, contribution)
      
      if (i == 1) {
        nrxn1_nlgn1_summary <- nrxn1_nlgn1_summary
      } else {
        nrxn1_nlgn1_summary <- bind_rows(nrxn1_nlgn1_summary, nrxn1_nlgn1_simple)
      }
    }
    
    # deal with Others sheet
    others_pairs_data <- tryCatch({
      read_excel(pairs_file, sheet = "Others")
    }, error = function(e) {
      cat("folder", folder, "Not found Pairs的Others sheet\n")
      return(NULL)
    })
    
    if (!is.null(others_pairs_data) && nrow(others_pairs_data) > 0) {
      others_pairs_simple <- others_pairs_data %>%
        select(name, contribution) %>%
        mutate(comb_id = comb_id) %>%
        select(comb_id, name, contribution)
      
      if (i == 1) {
        others_pairs_summary <- others_pairs_simple
      } else {
        others_pairs_summary <- bind_rows(others_pairs_summary, others_pairs_simple)
      }
    }
  }
}

# 3. Convert the data format and calculate the average
# 3.1 Convert attenuation.csv
attenuation_final <- attenuation_summary
write.csv(attenuation_final, "attenuation_summary.csv", row.names = FALSE)

# 3.2 Convert Counts_results.csv
# It's already in wide format. Save it directly.
counts_final <- counts_summary
write.csv(counts_final, "Counts_summary.csv", row.names = FALSE)


# 3.3 Convert Pathways.xlsx
# NRXN sheet
if (nrow(nrxn_summary) > 0) {
  nrxn_wide <- nrxn_summary %>%
    pivot_wider(
      id_cols = pathway,
      names_from = comb_id,
      values_from = ed
    )
  
  # Calculate the average value for each cycle
  avg_row <- data.frame(pathway = "average", t(colMeans(nrxn_wide[, -1], na.rm = TRUE)))
  names(avg_row) <- names(nrxn_wide)
  nrxn_with_avg <- bind_rows(nrxn_wide, avg_row)
  
  write.csv(nrxn_with_avg, "NRXN_summary.csv", row.names = FALSE)
}

# 3.4 Others sheet
if (nrow(others_summary) > 0) {
  others_wide <- others_summary %>%
    pivot_wider(
      id_cols = pathway,
      names_from = comb_id,
      values_from = ed
    )
  
  avg_row_others <- data.frame(pathway = "average", t(colMeans(others_wide[, -1], na.rm = TRUE)))
  names(avg_row_others) <- names(others_wide)
  others_with_avg <- bind_rows(others_wide, avg_row_others)
  
  write.csv(others_with_avg, "Others_pathways_summary.csv", row.names = FALSE)
}

# 3.5 Convert Pairs.xlsx - NRXN3_NLGN1 sheet
if (nrow(nrxn3_nlgn1_summary) > 0) {
  nrxn3_nlgn1_wide <- nrxn3_nlgn1_summary %>%
    pivot_wider(
      id_cols = name,
      names_from = comb_id,
      values_from = contribution
    )
  
  avg_row_nrxn3 <- data.frame(name = "average", t(colMeans(nrxn3_nlgn1_wide[, -1], na.rm = TRUE)))
  names(avg_row_nrxn3) <- names(nrxn3_nlgn1_wide)
  nrxn3_nlgn1_with_avg <- bind_rows(nrxn3_nlgn1_wide, avg_row_nrxn3)
  
  write.csv(nrxn3_nlgn1_with_avg, "NRXN3_NLGN1_summary.csv", row.names = FALSE)
}

# 3.6 Convert Pairs.xlsx - NRXN1_NLGN1 sheet
if (nrow(nrxn1_nlgn1_summary) > 0) {
  nrxn1_nlgn1_wide <- nrxn1_nlgn1_summary %>%
    pivot_wider(
      id_cols = name,
      names_from = comb_id,
      values_from = contribution
    )
  
  avg_row_nrxn1 <- data.frame(name = "average", t(colMeans(nrxn1_nlgn1_wide[, -1], na.rm = TRUE)))
  names(avg_row_nrxn1) <- names(nrxn1_nlgn1_wide)
  nrxn1_nlgn1_with_avg <- bind_rows(nrxn1_nlgn1_wide, avg_row_nrxn1)
  
  write.csv(nrxn1_nlgn1_with_avg, "NRXN1_NLGN1_summary.csv", row.names = FALSE)
}

# 3.7 Convert Pairs.xlsx - Others sheet
if (nrow(others_pairs_summary) > 0) {
  others_pairs_wide <- others_pairs_summary %>%
    pivot_wider(
      id_cols = name,
      names_from = comb_id,
      values_from = contribution
    )
  
  avg_row_others_pairs <- data.frame(name = "average", t(colMeans(others_pairs_wide[, -1], na.rm = TRUE)))
  names(avg_row_others_pairs) <- names(others_pairs_wide)
  others_pairs_with_avg <- bind_rows(others_pairs_wide, avg_row_others_pairs)
  
  write.csv(others_pairs_with_avg, "Others_pairs_summary.csv", row.names = FALSE)
}

cat("\nProcessing completed! The generated summary file has been generated：\n")
cat("1. attenuation_summary.csv\n")
cat("2. Counts_summary.csv\n")
cat("3. NRXN_summary.csv\n")
cat("4. Others_pathways_summary.csv\n")
cat("5. NRXN3_NLGN1_summary.csv\n")
cat("6. NRXN1_NLGN1_summary.csv\n")
cat("7. Others_pairs_summary.csv\n")




data_NRXN <- read.csv("NRXN_summary.csv", header=TRUE, row.names=1, check.names=FALSE)
rows_to_combine <- c("NRXN", "NRXN1", "NRXN2", "NRXN3", "NRXN4", "NRXN5", 
                     "NRXN6", "NRXN7", "NRXN8", "NRXN9")
if(!"row_names" %in% names(data_NRXN)) {
  data_NRXN$row_names <- rownames(data_NRXN)
}
data_NRXN_combine <- data_NRXN %>%
  filter(row_names %in% rows_to_combine) %>%  # 筛选需要的行
  select(-row_names) %>%  # 移除行名列
  t() %>%  # 转置
  as.vector() %>%  # 转换为向量
  # matrix(nrow = 1) %>%  # 转换为1行矩阵
  as.data.frame()

write.csv(data_NRXN_combine, "data_NRXN_combine.csv", row.names = FALSE)


data_NRXN_other <- read.csv("Others_pathways_summary.csv", header=TRUE, row.names=1, check.names=FALSE)

data_NRXN_other_combine <- pivot_longer(
  data = as.data.frame(data_NRXN_other), 
  cols = everything(),        
  names_to = "sample",        
  values_to = "value"       
)

data_NRXN_other_combine <- na.omit(data_NRXN_other_combine)
write.csv(data_NRXN_other_combine, "data_NRXN_other_combine.csv", row.names = FALSE)


data_pairs_others <- read.csv("Others_pairs_summary.csv", header=TRUE, row.names=1, check.names=FALSE)
data_pairs_others_without_avg <- data_pairs_others[-nrow(data_pairs_others), ]
data_pairs_others_without_avg <- pivot_longer(
  data = as.data.frame(data_pairs_others_without_avg),
  cols = everything(),        
  names_to = "sample",        
  values_to = "others"   
)

data_pairs_others_without_avg <- na.omit(data_pairs_others_without_avg)


write.csv(data_pairs_others_without_avg, "Others_pairs_combined.csv", row.names = FALSE)














