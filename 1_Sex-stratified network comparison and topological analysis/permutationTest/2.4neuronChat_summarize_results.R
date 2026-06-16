library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(openxlsx)

setwd("neuronChatF")
# 1. Get all result folders
iter_folders <- list.dirs(path = ".", recursive = FALSE, full.names = FALSE)
iter_folders <- iter_folders[grepl("^iter_\\d+_", iter_folders)]
cat("find", length(iter_folders), "result folders\n")

# Initialize the summary table
cell_comm_summary <- data.frame()

# 2. Traverse each folder
for (i in seq_along(iter_folders)) {
  folder <- iter_folders[i]
  cat("Dealing with folder:", folder, "...\n")
  
  # extract comb_id
  comb_id <- gsub("^iter_(\\d+)_.*", "\\1", folder)
  
  # file path
  result_file <- file.path(folder, "result.csv")
  
  # 2.2 Dealing with result.csv
  if (file.exists(result_file)) {
    result_data <- tryCatch({
      read.csv(result_file, header = TRUE, row.names = 1, check.names = FALSE)
    }, error = function(e) {
      cat("folder", folder, "Error reading result.csv\n")
      return(NULL)
    })
    
    if (!is.null(result_data)) {
      # Compute the emission value 
      # (sum of absolute values for each row): each row corresponds to the emission value of a cell
      # Compute the reception value 
      # (sum of absolute values for each column): each column corresponds to the reception value of a cell
      # Retrieve the cell names (row names and column names should be the same)
      cell_names <- rownames(result_data)
      outgoing_values <- apply(abs(result_data), 1, sum)
      incoming_values <- apply(abs(result_data), 2, sum)
      
      cell_comm_result <- data.frame(
        comb_id = comb_id,
        cell_type = rep(cell_names, 2),
        direction = c(rep("out", length(cell_names)), rep("in", length(cell_names))),
        value = c(outgoing_values, incoming_values)
      )
      
      cell_comm_wide <- data.frame(
        comb_id = comb_id,
        OPC1_out = outgoing_values["OPC1"],
        OPC2_out = outgoing_values["OPC2"],
        OPC3_out = outgoing_values["OPC3"],
        Oli1_out = outgoing_values["Oli1"],
        Oli2_out = outgoing_values["Oli2"],
        Oli3_out = outgoing_values["Oli3"],
        OPC1_in = incoming_values["OPC1"],
        OPC2_in = incoming_values["OPC2"],
        OPC3_in = incoming_values["OPC3"],
        Oli1_in = incoming_values["Oli1"],
        Oli2_in = incoming_values["Oli2"],
        Oli3_in = incoming_values["Oli3"]
      )
      
      # Append to the summary table
      if (i == 1) {
        cell_comm_summary <- cell_comm_wide
      } else {
        cell_comm_summary <- bind_rows(cell_comm_summary, cell_comm_wide)
      }
    }
  }
}

# 3. Convert the data format and calculate the average

# 3.1Save the cell communication results summary table.
if (nrow(cell_comm_summary) > 0) {
  write.csv(cell_comm_summary, "cell_communication_summary.csv", row.names = FALSE)
  cat("cell_communication_summary.csv\n")
  
  cell_comm_data <- cell_comm_summary
  
  cell_comm_combined <- cell_comm_data %>%
    mutate(
      OPC_out = OPC1_out + OPC2_out + OPC3_out,
      OPC_in = OPC1_in + OPC2_in + OPC3_in,
      OL_out = Oli1_out + Oli2_out + Oli3_out,
      OL_in = Oli1_in + Oli2_in + Oli3_in
    ) %>%
    select(comb_id, OPC_out, OPC_in, OL_out, OL_in)
  
  write.csv(cell_comm_combined, "cell_communication_combined_summary.csv", row.names = FALSE)
  cat("cell_communication_combined_summary.csv\n")
  
  
  cell_comm_stacked2 <- cell_comm_summary %>%
    pivot_longer(
      cols = -comb_id,
      names_to = "cell_direction",
      values_to = "value"
    ) %>%
    mutate(
      cell_type = substr(cell_direction, 1, 3),  
      index = substr(cell_direction, 4, 4),       
      direction = substr(cell_direction, 6, nchar(cell_direction))  #
    ) %>%
    
    pivot_wider(
      id_cols = c(comb_id, index),
      names_from = c(cell_type, direction),
      values_from = value
    ) %>%
    
    rename(
      OPC_out = OPC_out,
      OPC_in = OPC_in,
      OL_out = Oli_out,
      OL_in = Oli_in
    ) %>%
    arrange(comb_id, index)
  
  write.csv(cell_comm_stacked2, "cell_communication_stacked2.csv", row.names = FALSE)
  cat("cell_communication_stacked2.csv\n")
}

