#load libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(paletteer)
library(CellChat)
library(ggsignif)
library(stringr)
library(readxl)
library(ggrepel) 

#load SampleData 
Female_MDD <- readRDS("Female_MDD.rds")
Male_MDD <- readRDS("Male_MDD.rds")
Female_control <- readRDS("Female_control.rds")
Male_control <- readRDS("Male_control.rds")


#load dependencies code
source("source/scatter.R")#netAnalysis_signalingRole_scatter_V2
source("source/ChangesScatter.R")#netAnalysis_signalingChanges_scatter_V2
source("source/netAnalysis_contribution.R")#netAnalysis_contribution_v2


### attenuation analysis 

f_ca_mat <- Female_MDD@net$prob
f_ct_mat <- Female_control@net$prob
m_ca_mat <- Male_MDD@net$prob
m_ct_mat <- Male_control@net$prob
nn <- list(f_ct_mat,f_ca_mat,m_ct_mat,m_ca_mat)

compute_global_strength <- function(probability_matrix) { sum(abs(probability_matrix)) }

global_male_control <- compute_global_strength(m_ct_mat)
global_male_mdd <- compute_global_strength(m_ca_mat)
global_female_control <- compute_global_strength(f_ct_mat)
global_female_mdd <- compute_global_strength(f_ca_mat)

male_attenuation <- (1 - global_male_mdd / global_male_control) *100 # 0.04104659
female_attenuation <- (1 - global_female_mdd / global_female_control) *100 # 0.2303808


data <- data.frame(
  Gender = c("Male", "Female"),
  Attenuation = c(0.04104659*100, 0.2303808*100),
  Group = "MDD vs Control"
)

# F1a 
ggplot(data, aes(x = Gender, y = Attenuation, fill = Gender)) +
  geom_bar(stat = "identity", width = 0.3, alpha = 0.8, position = position_dodge(width = 0.05) ) +
  geom_text(aes(label = sprintf("%.2f%%", Attenuation)), vjust = -0.5, size = 5, color = "black") +  
  scale_fill_manual(values = c("Female" = "#D62728","Male" = "#4682B4")) +  
  labs(
    title = "Signal Attenuation in MDD vs. Control\n (1 - MDD/Control)",
    x = "Gender",
    y = "Attenuation (%)"
  ) +
  theme_minimal(base_size = 12) +
  coord_cartesian(ylim = c(0, 30)) +
  theme(
    plot.title = element_text(size = 13, hjust = 0.5),
    plot.margin = margin(5, 5, 5, 5, "mm"),
    aspect.ratio = 0.4,
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.y = element_blank(),
    legend.position = "none",  
    axis.text.x = element_text(face = "bold", size = 10)
  )


#F1c

cellchat_list <- readRDS("cellchat_list.rds")

female_case_list <- c("F1", "F11", "F12", "F14", "F15", "F16", "F17", "F18", "F19", 
                      "F2", "F20", "F25", "F27", "F28", "F3", "F4", "F5", "F6","F8","F9")

female_control_list <- c("F10", "F13", "F21", "F22", "F23", "F24", "F26", "F29", 
                         "F30", "F31", "F32", "F33", "F34", "F35", "F36", "F37", 
                         "F38", "F7")

male_case_list <- c("M1", "M10", "M11", "M14", "M17", "M18", "M23", "M26", "M28", "M30", "M32", "M33", "M34", "M4", "M5", "M6", "M8")
male_control_list <- c("M12", "M13", "M15", "M16", "M19", "M2", "M20", "M21", "M22", "M24", "M27", "M29", "M3", "M31", "M7", "M9")


sample_info <- list()
for (sample in female_case_list) {
  sample_info[[sample]] <- c(gender = "Female", status = "Case")
}
for (sample in female_control_list) {
  sample_info[[sample]] <- c(gender = "Female", status = "Control")
}
for (sample in male_case_list) {
  sample_info[[sample]] <- c(gender = "Male", status = "Case")
}
for (sample in male_control_list) {
  sample_info[[sample]] <- c(gender = "Male", status = "Control")
}


cell_strength_df <- data.frame(
  Sample = character(),
  CellType = character(),
  Outgoing_Strength = numeric(),
  Incoming_Strength = numeric(),
  Gender = character(),
  Status = character(),
  Group = character(),
  stringsAsFactors = FALSE
)

# Calculate the outgoing and incoming signaling strength for each cell type in each sample.
# Ast1         Ast2         End1     ExN1_L24    ExN10_L46    ExN11_L56    ExN12_L56    ExN13_L56        ExN14 
# 166           12          264           25           39           17           23           11           31 
# ExN15_L56    ExN16_L56        ExN17        ExN18    ExN19_L56     ExN2_L23    ExN20_L56     ExN3_L46     ExN4_L35 
# 13           15           80           10            4          529            1          127          443 
# ExN5         ExN6         ExN7     ExN8_L24     ExN9_L23      InN1_PV InN10_ADARB2     InN2_SST     InN3_VIP 
# 21           61           76           81           44          155           14           51           32 
# InN4_VIP     InN5_SST   InN6_LAMP5     InN7_Mix     InN8_Mix      InN9_PV         Mic1          Mix         Oli1 
# 147          133           48           18           42           23          252           99           48 
# Oli2         Oli3         OPC1         OPC2         OPC3 
# 30          406          146          125            1 
for (sample in names(sample_info)) {
  if (sample %in% names(cellchat_list)) {
    weight_matrix <- cellchat_list[[sample]]@net$weight
    
    cell_types <- rownames(weight_matrix)
    
    outgoing_strength <- colSums(weight_matrix, na.rm = TRUE)
    
    incoming_strength <- rowSums(weight_matrix, na.rm = TRUE)
    
    for (i in 1:length(cell_types)) {
      cell_strength_df <- rbind(cell_strength_df, data.frame(
        Sample = sample,
        CellType = cell_types[i],
        Outgoing_Strength = outgoing_strength[i],
        Incoming_Strength = incoming_strength[i],
        Gender = sample_info[[sample]][1],
        Status = sample_info[[sample]][2],
        stringsAsFactors = FALSE
      ))
    }
  } else {
    warning(paste("sample", sample, " not in cellchat_list"))
  }
}

cell_strength_df$Total_Strength <- cell_strength_df$Outgoing_Strength + cell_strength_df$Incoming_Strength

cell_strength_df$Group <- paste(cell_strength_df$Gender, cell_strength_df$Status, sep = "_")
cell_strength_df$Group <- gsub("Case", "MDD", cell_strength_df$Group)
cell_strength_df$Group <- factor(cell_strength_df$Group, 
                                 levels = c("Female_MDD", "Female_Control", 
                                            "Male_MDD", "Male_Control"))


cell_strength_long <- pivot_longer(
  data = cell_strength_df,
  cols = c("Outgoing_Strength", "Incoming_Strength"), 
  names_to = "Strength_Type",  
  values_to = "Strength"       
)


colors <- c(
  "Female_MDD"     = "#D62728",  
  "Female_Control" = "#F7B7B0",  
  "Male_MDD"       = "#4682B4",  
  "Male_Control"   = "#A3C4DC"  
)


ggplot(cell_strength_df, aes(x = Group, y = Total_Strength, fill = Group)) +
  geom_violin(alpha = 0.7, linewidth = 1.2) +
  scale_fill_manual(values = colors) +
  geom_boxplot(width = 0.1, outlier.alpha = 0.5, linewidth = 0.8) +
  theme_classic(base_size = 20) +
  theme(
    axis.text = element_text(color = 'black'),
    legend.position = 'none'
  ) +
  labs(
    x = "",
    y = ""
  )


#F1d
gender <- c("Female","Male")
direc <- c("Transmission","Reception")

get_hub_strength <- function(prob_matrix, dim =1,top=0.3) {
  strength <- apply(prob_matrix, dim, sum)
  sorted_idx <- order(strength, decreasing = TRUE)
  n_hubs <- round(top * length(strength))
  hub_global_idx <- sorted_idx[1:n_hubs]
  hub_strength <- strength[hub_global_idx]
  return(hub_strength)
}

n_iter <- length(gender) * length(direc)

gender_labels <- vector("character", length = n_iter)
direc_labels <- vector("character", length = n_iter)
counter <- 1
rank_df <- list()

for (i in 1:length(gender)) {
  for (j in 1:length(direc)) {
    hub_ctrl <- get_hub_strength(nn[[2*i-1]],j)
    hub_mdd <- get_hub_strength(nn[[2*i]],j)
    
    common_cells <- intersect(names(hub_ctrl), names(hub_mdd))
    rank_ctrl <- match(common_cells, names(hub_ctrl))
    rank_mdd <- match(common_cells, names(hub_mdd))
    common_ctrl <- hub_ctrl[common_cells]
    common_mdd <- hub_mdd[common_cells]
    
    gender_labels[counter] <- gender[i]
    direc_labels[counter] <- direc[j]
    a <- data.frame(
      Gender = gender[i],
      Direction = direc[j],
      CellType = common_cells,
      Ctrl_Rank = rank_ctrl,
      MDD_Rank = rank_mdd,
      Rank_Diff = rank_ctrl - rank_mdd
    )
    rank_df[[counter]] <-a
    counter <- counter + 1
  }
}

rank_df <- do.call(rbind, rank_df)

rank_df <- rank_df %>%
  mutate(
    Gender = factor(Gender, levels= gender),
    Direction = factor(Direction, levels = direc)
  )

rank_df <- rank_df %>%
  mutate(
    Rank_Change = case_when(
      Rank_Diff > 0 ~ "Improved",
      Rank_Diff < 0 ~ "Declined",
      TRUE ~ "Unchanged"),
    Rank_Change_Factor = factor(Rank_Change, levels = c("Improved", "Unchanged", "Declined")) )

rank_df <- rank_df %>% 
  mutate(
    Group_Label = interaction(rank_df$Gender, rank_df$Direction,sep = "_")
  )
rank_df <- rank_df %>% 
  mutate(
    Group_Label = factor(Group_Label, levels = levels(Group_Label))
  )

ggplot(rank_df, aes(x = Ctrl_Rank, y = MDD_Rank, color = Rank_Change_Factor)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_segment(aes(xend = Ctrl_Rank, yend = Ctrl_Rank), color = "grey50", linetype = "dashed") +
  geom_abline(slope = 1, intercept = 0, color = "#228B22", linetype = "solid") +
  ggrepel::geom_text_repel(
    aes(label = CellType),
    size = 6,
    max.overlaps = 20,
    box.padding = 0.5
  ) +
  facet_wrap(~ Group_Label, ncol = 2, scales = "free") +
  scale_color_manual(
    values = c("Improved" = "#6E1633", "Unchanged" = "#4D3015", "Declined" = "#0C6736"),
    name = "Rank Change"
  ) +
  labs(
    title = "Cell Type Rank Changes in MDD",
    subtitle = "Comparison Across Gender and Direction",
    x = "Control Group Rank",
    y = "MDD Group Rank"
  ) +
  theme_bw(base_size = 12) +coord_cartesian(xlim = c(0, 12.5)) +coord_cartesian(ylim = c(0, 12.5)) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom"
  )



#F1e

object.list.all <- list(
  Female_control = Female_control,
  Female_MDD = Female_MDD,
  Male_control = Male_control,
  Male_MDD = Male_MDD
)
cellchat <- mergeCellChat(object.list.all, add.names = names(object.list.all))

netAnalysis_signalingRole_scatter_V2(cellchat, data=c("Male_control","Male_MDD","Female_control","Female_MDD"), opc_flag = F, opc_flag2 = F, gender_flag = F, gender_flag2 = T, F_flag = F, M_flag= F)


#F1g



##################################################################################
#F2a
df <- expand.grid(
  cell_type = c("ExN","InN","OPC","Oli"),                 
  sex = c("Females", "Males"),               
  layer = c("Euclidean Distance", "Horizontal Distance", "Vertical Distance") 
)

df <- df %>%
  mutate(
    layer_label = case_when(
      layer == "Euclidean Distance" ~ "Euclidean Distance",
      layer == "Horizontal Distance" ~ "Horizontal Distance",
      layer == "Vertical Distance" ~ "Vertical Distance"
    )
  )

#df$Count <- c(8,1,1,0,7,1,1,1,8,2,0,0,9,0,1,0,6,1,2,1,4,2,2,2)  # male + female
df$Count <- c(7,1,1,1,8,1,1,0,9,0,1,0,8,2,0,0,4,2,2,2,6,1,2,1)  # female + male
df$Total <- rep(c(21,10,3,3,20,10,3,3), 3)                          # No. of superclass (fixed)
df$Ratio <- df$Count / df$Total

ggplot(df,aes(x=sex)) +
  geom_col(
    aes(y = 1,
        group = cell_type),
    fill = "grey80",
    width = 0.7,
    position = position_dodge(width = 0.8),
    alpha = 0.4
  ) +
  geom_col(
    aes(y = Ratio,
        fill = cell_type),
    width = 0.6,
    position = position_dodge(width = 0.8),
    color = "black",
    linewidth = 0.3
  ) +
  facet_wrap(~ layer_label, nrow = 1 ) +
  scale_fill_manual(
    values = c("#DC9CC3", "#A7D5B6", "#E4A71D", "#D6C576"),
    name = "Cell Types"
  ) +
  labs(
    title = "The Most Interactive Cells Between MDD Patients and Controls",
    x = "Gender",
    y = "Ratio"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0),
    panel.spacing = unit(1.8, "lines"),
    #strip.background = element_rect(fill = "grey95"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.4, size = 16)
  )

#######################f2c

object.list_F <- list(
  Female_control = Female_control,
  Female_MDD = Female_MDD
)
cellchat_F <- mergeCellChat(object.list_F, add.names = names(object.list_F))
netAnalysis_signalingChanges_scatter_V2(cellchat_F,idents.use.multi =c("ExN7", "OPC1", "ExN17", "InN3_VIP", "ExN8_L24", "ExN14", "ExN13_L56", "Oli1", "ExN18", "ExN11_L56"))

###########################f2d

object.list_M <- list(
  Male_control = Male_control,
  Male_MDD = Male_MDD
)
cellchat_M <- mergeCellChat(object.list_M, add.names = names(object.list_M))
netAnalysis_signalingChanges_scatter_V2(cellchat_M,idents.use.multi =  c("ExN12_L56", "ExN13_L56", "ExN15_L56", "ExN16_L56", "ExN18", "ExN20_L56", "ExN6", "ExN7", "InN7_Mix", "OPC3"))


###########################f2f
object.list.all <- list(
  Female_control = Female_control,
  Female_MDD = Female_MDD,
  Male_control = Male_control,
  Male_MDD = Male_MDD
)
cellchat <- mergeCellChat(object.list.all, add.names = names(object.list.all))

netAnalysis_contribution_v2(cellchat, dataset = c("Female_control","Female_MDD","Male_control","Male_MDD"), signaling = "NRXN")



###########################f2f
#The siggene are derived from other articles.
file <- read_excel("../data/CSF.xlsx")
head(file)
colnames(file) <- c("Gene", "NegLogP", "QValue", "FoldChange", "TestStatistic")
file$FoldChange <- as.numeric(as.character(file$FoldChange))
file$NegLogP <- as.numeric(as.character(file$NegLogP))
file$Significant <- ifelse(file$NegLogP > -log10(0.05) & abs(file$FoldChange) > 0.5, "Significant", "Not Significant")
file$Score <- file$NegLogP * abs(file$FoldChange)

top10_genes <- file[order(-file$Score), ]  
top10_genes <- top10_genes[1:12, ]     
target_genes <- c("NRXN1", "NRXN2", "NRXN3", "NLGN1")
file$Highlight <- ifelse(file$Gene %in% target_genes, "Highlighted", NA)

ggplot(file, aes(x = FoldChange, y = NegLogP, color = Significant)) +
  geom_point(alpha = 0.8, size = 2) + 
  geom_point(
    data = subset(file, !is.na(Highlight)), 
    aes(color = "Highlighted"), 
    size = 3, 
    shape = 16 
  ) +  
  scale_color_manual(
    values = c("Significant" = "#00A087", "Not Significant" = "#3C5488", "Highlighted" = "#F39B7F"),
    name = "Significance"
  ) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "#B22222") +  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#B22222") + 
  geom_text_repel(
    data = top10_genes,
    aes(label = Gene),
    size = 4,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = "grey50"
  ) +  
  geom_text_repel(
    data = subset(file, !is.na(Highlight)), 
    aes(label = Gene, color = "Highlighted"),  
    size = 5,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = "#B22222"
  ) +  
  labs(
    title = "Volcano Plot",
    x = "Log2(Fold Change)",
    y = "-Log10(P-value)",
    color = "Significance"
  ) +
  theme_minimal()


#S2a/b
netVisual_heatmap(cellchat_F, measure = "weight",font.size = 12,font.size.title = 16)
netVisual_heatmap(cellchat_M, measure = "weight",font.size = 12,font.size.title = 16)

