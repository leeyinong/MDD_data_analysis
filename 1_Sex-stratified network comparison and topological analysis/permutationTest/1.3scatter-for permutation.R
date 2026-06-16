netAnalysis_signalingRole_scatter_V2 <- function (object, signaling = NULL, color.use = NULL, slot.name = "netP", 
                                                  group = NULL, weight.MinMax = NULL, dot.size = c(1, 5), 
                                                  point.shape = c(21, 22, 24, 23, 25, 8, 3), label.size = 3, 
                                                  dot.alpha = 0.9, x.measure = "outdeg", y.measure = "indeg", 
                                                  dataset = c(NULL,NULL), legend.shape = c(16,15,17,18),
                                                  xlabel = "Outgoing interaction strength", 
                                                  ylabel = "Incoming interaction strength", 
                                                  title = NULL, font.size = 10, font.size.title = 10, opc_flag = F, gender_flag = F, F_flag = F, M_flag = F,
                                                  gender_flag2 = F, opc_flag2 = F, do.label = TRUE, show.legend = TRUE, show.axes = TRUE,external_objects = NULL) 
{
  #==============================================================================================
  for (i in 1:length(dataset)) {
    if (length(slot(object, slot.name)[[dataset[i]]]$centr) == 0) {
      stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
    }
    
    if (sum(c(x.measure, y.measure) %in% names(slot(object, slot.name)[[dataset[i]]]$centr[[1]])) != 2) {
      
      stop(paste0("`x.measure, y.measure` should be one of ", 
                  paste(names(slot(object, slot.name)[[dataset[i]]]$centr), 
                        collapse = ", "), 
                  "\n", "`outdeg_unweighted` is only supported for version >= 1.1.2"))
    }
  }
  
  cent <- list()
  outgoing <- list()
  incoming <- list()
  outgoing.cells <- list()
  incoming.cells <- list()
  num.link <- list()
  stre.link <- list()
  slope <- list()
  ed <- list()
  data <-""
  gg<-NULL
  
  for (i in 1:length(dataset)) {
    
    cent[[i]] <- methods::slot(object, slot.name)[[dataset[i]]]$centr
    outgoing[[i]] <- matrix(0, nrow = nlevels(object@idents[[dataset[i]]]), ncol = length(cent[[i]]))
    incoming[[i]] <- matrix(0, nrow = nlevels(object@idents[[dataset[i]]]), ncol = length(cent[[i]]))
    name <- paste0(levels(object@idents[[dataset[i]]]),"--",dataset[i])
    dimnames(outgoing[[i]]) <- list(name, names(cent[[i]]))
    dimnames(incoming[[i]]) <- dimnames(outgoing[[i]])
    
    for (j in 1:length(cent[[i]])) {
      outgoing[[i]][, j] <- cent[[i]][[j]][[x.measure]]    # outgoing[[i]]: cell clusters * pathways
      incoming[[i]][, j] <- cent[[i]][[j]][[y.measure]]
    }
    
    if (is.null(signaling)) {
      message("Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways")
    } else {
      message("Signaling role analysis on the cell-cell communication network from user's input")
      signalings <- signaling[signaling %in% object@netP[[dataset[i]]]$pathways]
      
      if (length(signaling) == 0) {
        stop("There is no significant communication for the input signaling. All the significant signaling are shown in `object@netP$pathways`")
      }
      outgoing[[i]] <- outgoing[[i]][, signalings, drop = FALSE]
      incoming[[i]] <- incoming[[i]][, signalings, drop = FALSE]
    }
    
    outgoing.cells[[i]] <- rowSums(outgoing[[i]])
    incoming.cells[[i]] <- rowSums(incoming[[i]])
    
    #######################################################
    
    # tmp <- paste0(dataset[i],".rds")
    # data <- readRDS(tmp)
    
    ######################################################
    obj <- aggregateNet(external_objects[[dataset[i]]], signaling = signaling, remove.isolate = FALSE)
    num.link[[i]] <- obj@net$count
    stre.link[[i]] <- obj@net$weight
    
    num.link[[i]] <- rowSums(num.link[[i]]) + colSums(num.link[[i]]) - diag(num.link[[i]])
    stre.link[[i]] <- rowSums(stre.link[[i]]) + colSums(stre.link[[i]]) - diag(stre.link[[i]])
    
  }
  outgoing.cells <- unlist(outgoing.cells)
  incoming.cells <- unlist(incoming.cells)
  num.link <- unlist(num.link)
  stre.link <- unlist(stre.link)
  slope <- ifelse(outgoing.cells != 0, incoming.cells/outgoing.cells, 0)
  ed <- sqrt((outgoing.cells)^2+(incoming.cells)^2)
  
  group <- sub(".*--", "", names(incoming.cells))
  labels <- sub("--.*", "", names(incoming.cells))
  
  # the levels of clusters in the merged object are fixed by the levels of clusters in the female_case dataset
  df <- data.frame(x = outgoing.cells, y = incoming.cells, slope = slope, ed = ed,
                   labels = factor(labels, levels = levels(object@idents[[2]])), Interactions = num.link, Strength=stre.link, Group = factor(group, levels = unique(group)))

  # write.csv(df,"signalingRole.csv")
  
  df <- df %>% mutate(Condition = factor(ifelse(str_detect(Group, "MDD"), "MDD", "Control"),levels = c("MDD", "Control")))
  df <- df %>% mutate(Female = ifelse(str_detect(Group, "Female"), ifelse(str_detect(Group, "MDD"), "MDD", "Control"), "Others"))
  df <- df %>% mutate(Male = ifelse(str_detect(Group, "Male"), ifelse(str_detect(Group, "MDD"), "MDD", "Control"), "Others"))
 
  df$Cells_Types <- case_when(
    str_detect(labels, "Ast") ~ "Ast" ,
    str_detect(labels, "End") ~ "End" ,
    str_detect(labels, "ExN") ~ "ExN" ,
    str_detect(labels, "InN") ~ "InN" ,
    str_detect(labels, "Mic") ~ "Mic" ,
    str_detect(labels, "Mix") ~ "Mix" ,
    str_detect(labels, "OPC") ~ "OPC" ,
    str_detect(labels, "Oli") ~ "OL" 
  )
  
  # write.csv(df,"D:/Project/niu/OPC-MDD/permutation/test.csv")
  write.csv(df,"test.csv")
  
  df_diff <- df %>%
    pivot_wider(
      id_cols = labels,  
      names_from = Condition,
      values_from = c(x, y),
      names_sep = "_"
    ) %>%
    mutate(
      x_diff = abs(x_MDD - x_Control),
      y_diff = abs(y_MDD - y_Control),
      ed_diff = sqrt((x_diff)^2+(y_diff)^2)
    )

  df_x <- df_diff %>% select(labels,x_diff) %>% arrange(desc(x_diff))
  colnames(df_x) <- c("cell_type","sth")
  top_10_x <- head(df_x, 10)
  #write.csv(df_x,"D:/Project/niu/OPC-MDD/permutation/test_x.csv", row.names = FALSE)
  df_y <- df_diff %>% select(labels,y_diff) %>% arrange(desc(y_diff))
  colnames(df_y) <- c("cell_type","sth")
  top_10_y <- head(df_y, 10)
  #write.csv(df_y,"D:/Project/niu/OPC-MDD/permutation/test_y.csv", row.names = FALSE)
  df_ed <- df_diff %>% select(labels,ed_diff) %>% arrange(desc(ed_diff))
  colnames(df_ed) <- c("cell_type","sth")
  top_10_ed <- head(df_ed, 10)
  #write.csv(df_ed,"D:/Project/niu/OPC-MDD/permutation/test_ed.csv", row.names = FALSE)
  
  wb <- createWorkbook()
  
  addWorksheet(wb, "Differences")
  addWorksheet(wb, "Euclidean_Differences")
  addWorksheet(wb, "Horizontal_Differences")
  addWorksheet(wb, "Vertical_Differences")
  addWorksheet(wb, "Euclidean_Top10")
  addWorksheet(wb, "Horizontal_Top10")
  addWorksheet(wb, "Vertical_Top10")
  
  writeData(wb, sheet = "Differences", df_diff)
  writeData(wb, sheet = "Euclidean_Differences", df_ed)
  writeData(wb, sheet = "Horizontal_Differences", df_x)
  writeData(wb, sheet = "Vertical_Differences", df_y)
  writeData(wb, sheet = "Euclidean_Top10", top_10_ed)
  writeData(wb, sheet = "Horizontal_Top10", top_10_x)
  writeData(wb, sheet = "Vertical_Top10", top_10_y)
  
  # saveWorkbook(wb, "D:/Project/niu/OPC-MDD/permutation/Differences.xlsx", overwrite = TRUE)
  saveWorkbook(wb, "Differences.xlsx", overwrite = TRUE)
  
  }

