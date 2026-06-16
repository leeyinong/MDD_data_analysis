#load libraries
library(dplyr)
library(NeuronChat)
library(CellChat)
library(Seurat)
library(ggplot2)
library(ggalluvial)
library(ComplexHeatmap)
library(circlize)

#load SampleData 
fca_x <- readRDS("female_case_neuronchat.rds")
fct_x <- readRDS("female_control_neuronchat.rds")
mca_x <- readRDS("male_case_neuronchat.rds")
mct_x <- readRDS("male_control_neuronchat.rds")

#load dependencies code
source("source/heatmap_aggregated_v2.R")#heatmap_aggregated_v2
source("source/net_aggregation.R")#net_aggregation_v2



# group info
cluster_to_group <- data.frame(
  Cluster = c("Ast1", "Ast2", "End1", "ExN1_L24","ExN10_L46","ExN11_L56","ExN12_L56","ExN13_L56","ExN14","ExN15_L56","ExN16_L56","ExN17","ExN18","ExN19_L56","ExN2_L23","ExN20_L56","ExN3_L46","ExN4_L35","ExN5","ExN6","ExN7", "ExN8_L24","ExN9_L23","InN1_PV","InN10_ADARB2","InN2_SST","InN3_VIP","InN4_VIP","InN5_SST","InN6_LAMP5","InN7_Mix","InN8_Mix","InN9_PV","Mic1","Mix","Oli1","Oli2","Oli3","OPC1","OPC2","OPC3"),
  Group = c("Non-Neuronal", "Non-Neuronal", "Non-Neuronal", "Glutamatergic", "Glutamatergic", "Glutamatergic", "Glutamatergic","Glutamatergic", "Glutamatergic", "Glutamatergic", "Glutamatergic","Glutamatergic", "Glutamatergic", "Glutamatergic", "Glutamatergic","Glutamatergic", "Glutamatergic", "Glutamatergic", "Glutamatergic","Glutamatergic", "Glutamatergic", "Glutamatergic", "Glutamatergic","GABAergic","GABAergic","GABAergic","GABAergic","GABAergic","GABAergic","GABAergic","GABAergic","GABAergic","GABAergic","Non-Neuronal", "Non-Neuronal", "Non-Neuronal","Non-Neuronal", "Non-Neuronal", "Non-Neuronal","Non-Neuronal", "Non-Neuronal")
)
group <- setNames(cluster_to_group$Group, cluster_to_group$Cluster)

#F3a
data_list <- c("Female_Control","Female_MDD")
nc <- mergeNeuronChat(list(fct_x,fca_x), add.names = data_list)
heatmap_aggregated_v2(nc, dataset=data_list,method='weight',group = group)

#F3b
data_list <- c("Male_Control","Male_MDD")
nc <- mergeNeuronChat(list(mct_x,mca_x), add.names = data_list)
#Need to switch the matrix dimensions in the heatmap_aggregated_v2.R.(net_agg[36:41,36:41])
heatmap_aggregated_v2(nc, dataset=data_list,method='weight',group = group)

#F3d

mca_n1 <- mca_x@net[["NRXN1_NLGN1"]]
mct_n1 <- mct_x@net[["NRXN1_NLGN1"]]
mca_n3 <- mca_x@net[["NRXN3_NLGN1"]]
mct_n3 <- mct_x@net[["NRXN3_NLGN1"]]
fca_n1 <- fca_x@net[["NRXN1_NLGN1"]]
fct_n1 <- fct_x@net[["NRXN1_NLGN1"]]
fca_n3 <- fca_x@net[["NRXN3_NLGN1"]]
fct_n3 <- fct_x@net[["NRXN3_NLGN1"]]





nn <- list(mct_n3,mca_n3,fct_n3,fca_n3) # for NRXN3_NLGN1
val <- data.frame(Data=numeric(), Cluster= numeric(), Value=numeric(),stringsAsFactors = F)
for (i in 1: (length(nn)-2)){   # calculating male datasets
  print(i)
  sender <- nn[[i]][35:40,]
  receiver <- nn[[i]][,35:40]
  
  print(rowSums(sender))
  print(colSums(receiver))
  
  tmp <- data.frame(Data=i,Cluster=names(rowSums(sender)), Value=rowSums(sender))
  val <- rbind(val,tmp)
  tmp <- data.frame(Data=i,Cluster=names(colSums(receiver)), Value=colSums(receiver))
  val <- rbind(val,tmp)
}
for (i in (length(nn)-1): length(nn)){   # calculating female datasets
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

n_datasets <- 4
n_sexes <- 2
n_celltypes <- 6

data <- expand.grid(
  CellType = c("OPC1","OPC2","OPC3","Oli1","Oli2","Oli3"),
  Direction = c("Outgoing", "Incoming"),
  Dataset = c("Control","MDD"),
  Gender = c("Male","Female")
) %>%
  mutate(
    CellType = factor(CellType, levels = unique(CellType)),
    Gender = factor(Gender, levels = c("Male", "Female")),
    Dataset = factor(Dataset, levels = c("Control", "MDD")),
    Value = val$Value)

data <- data %>% mutate(
  Value = ifelse(Direction == "Outgoing", -1 * abs(Value), abs(Value)))

data <- data %>% mutate(
  GroupID = interaction(
    CellType, Direction, # CellType
    Dataset,Gender,
    sep = "_",
    lex.order = FALSE 
  )
) 

data <- data %>% mutate(
  Show = interaction(Dataset,Gender,sep = "_")
)

head(data)
# coordinates
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

ggplot(data, aes(x = x_final, y = Value)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.8) +
  geom_point(
    data = ~ filter(.x, Value != 0),  
    aes(shape = Dataset, color= Gender, fill = Gender),
    size = 3,
    stroke = 0.5
  ) +
  scale_x_continuous(
    breaks = unique(data$x_base),
    labels = levels(factor(data$CellType))
  ) +
  annotate(
    "text",
    x = min(data$x_final) - 0.4,                        
    y = 5,                
    label = "Incoming", 
    angle = 90,                      
    #hjust = c(1, -1),                 
    #vjust = c(1.5, -0.5),            
    size = 4,                        
    color = "black"
  )+
  annotate(
    "text",
    x = min(data$x_final) - 0.4,                        
    y = -5,                
    label = "Outgoing", 
    angle = 90,                       
    #hjust = c(1, -1),                
    #vjust = c(1.5, -0.5),           
    size = 4,                       
    color = "black"
  ) + coord_cartesian(xlim = c(0.5, 6.5), clip = "off")+
  scale_y_continuous(
    labels = abs,
    limits = c(-max(abs(data$Value)) * 1.1, max(abs(data$Value)) * 1.1),
    expand = expansion(mult = 0.1)
  ) +
  scale_color_manual(values = c("Female" = "#d62728", "Male" = "steelblue")) +
  scale_fill_manual(values = c("Female" = "#d62728", "Male" = "steelblue")) +
  scale_shape_manual(values = c("Control" = 22, "MDD" = 21)) + 
  guides(
    color = guide_legend(override.aes = list(fill = NA)), 
    shape = guide_legend(override.aes = list(size = 3, shape = c(22, 21)))  # Control=22, MDD=21
  )+ 
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.box = "vertical",
    legend.spacing.y = unit(0.5, "cm")
  ) +
  labs(
    x = "Cell Types",
    y = "Signal Strength",
    color = "Gender",
    fill = "Gender",
    shape = "Condition"
  )


#f3e

nn <- list(mct_n1,mca_n1,fct_n1,fca_n1) # for NRXN1_NLGN1
val <- data.frame(Data=numeric(), Cluster= numeric(), Value=numeric(),stringsAsFactors = F)
for (i in 1: (length(nn)-2)){   # calculating male datasets
  print(i)
  sender <- nn[[i]][35:40,]
  receiver <- nn[[i]][,35:40]
  
  print(rowSums(sender))
  print(colSums(receiver))
  
  tmp <- data.frame(Data=i,Cluster=names(rowSums(sender)), Value=rowSums(sender))
  val <- rbind(val,tmp)
  tmp <- data.frame(Data=i,Cluster=names(colSums(receiver)), Value=colSums(receiver))
  val <- rbind(val,tmp)
}
for (i in (length(nn)-1): length(nn)){   # calculating female datasets
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

n_datasets <- 4
n_sexes <- 2
n_celltypes <- 6

data <- expand.grid(
  CellType = c("OPC1","OPC2","OPC3","Oli1","Oli2","Oli3"),
  Direction = c("Outgoing", "Incoming"),
  Dataset = c("Control","MDD"),
  Gender = c("Male","Female")
) %>%
  mutate(
    CellType = factor(CellType, levels = unique(CellType)),
    Gender = factor(Gender, levels = c("Male", "Female")),
    Dataset = factor(Dataset, levels = c("Control", "MDD")),
    Value = val$Value)

data <- data %>% mutate(
  Value = ifelse(Direction == "Outgoing", -1 * abs(Value), abs(Value)))

data <- data %>% mutate(
  GroupID = interaction(
    CellType, Direction, # CellType
    Dataset,Gender,
    sep = "_",
    lex.order = FALSE 
  )
) 

data <- data %>% mutate(
  Show = interaction(Dataset,Gender,sep = "_")
)

head(data)

# coordinates
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
ggplot(data, aes(x = x_final, y = Value)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.8) +
  geom_point(
    data = ~ filter(.x, Value != 0),  
    aes(shape = Dataset, color= Gender, fill = Gender),
    size = 3,
    stroke = 0.5
  ) +
  scale_x_continuous(
    breaks = unique(data$x_base),
    labels = levels(factor(data$CellType))
  ) +
  annotate(
    "text",
    x = min(data$x_final) - 0.4,                        
    y = 5,                
    label = "Incoming", 
    angle = 90,                      
    #hjust = c(1, -1),                 
    #vjust = c(1.5, -0.5),            
    size = 4,                        
    color = "black"
  )+
  annotate(
    "text",
    x = min(data$x_final) - 0.4,                        
    y = -5,                
    label = "Outgoing", 
    angle = 90,                       
    #hjust = c(1, -1),                
    #vjust = c(1.5, -0.5),           
    size = 4,                       
    color = "black"
  ) + coord_cartesian(xlim = c(0.5, 6.5), clip = "off")+
  scale_y_continuous(
    labels = abs,
    limits = c(-max(abs(data$Value)) * 1.1, max(abs(data$Value)) * 1.),
    expand = expansion(mult = 0.1)
  ) +
  scale_color_manual(values = c("Female" = "#d62728", "Male" = "steelblue")) +
  scale_fill_manual(values = c("Female" = "#d62728", "Male" = "steelblue")) +
  scale_shape_manual(values = c("Control" = 22, "MDD" = 21)) + 
  guides(
    color = guide_legend(override.aes = list(fill = NA)), 
    shape = guide_legend(override.aes = list(size = 3, shape = c(22, 21)))  # Control=22, MDD=21
  )+ 
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.box = "vertical",
    legend.spacing.y = unit(0.5, "cm")
  ) +
  labs(
    x = "Cell Types",
    y = "Signal Strength",
    color = "Gender",
    fill = "Gender",
    shape = "Condition"
  )


#S4a/b/c/d
heatmap_aggregated(fct_x, method='weight',group = group)
heatmap_aggregated(fca_x, method='weight',group = group)
heatmap_aggregated(mct_x, method='weight',group = group)
heatmap_aggregated(mca_x, method='weight',group = group)

#S5
heatmap_single(fct_x,interaction_name='NRXN3_NLGN1',group=group)
heatmap_single(fca_x,interaction_name='NRXN3_NLGN1',group=group)
heatmap_single(mct_x,interaction_name='NRXN3_NLGN1',group=group)
heatmap_single(mca_x,interaction_name='NRXN3_NLGN1',group=group)

heatmap_single(fct_x,interaction_name='NRXN1_NLGN1',group=group)
heatmap_single(fca_x,interaction_name='NRXN1_NLGN1',group=group)
heatmap_single(mct_x,interaction_name='NRXN1_NLGN1',group=group)
heatmap_single(mca_x,interaction_name='NRXN1_NLGN1',group=group)
