netAnalysis_signalingRole_scatter_V2 <- function (object, signaling = NULL, color.use = NULL, slot.name = "netP", 
                                                  group = NULL, weight.MinMax = NULL, dot.size = c(1, 5), 
                                                  point.shape = c(21, 22, 24, 23, 25, 8, 3), label.size = 3, 
                                                  dot.alpha = 0.9, x.measure = "outdeg", y.measure = "indeg", 
                                                  dataset = c(NULL,NULL), legend.shape = c(16,15,17,18),
                                                  xlabel = "Outgoing interaction strength", 
                                                  ylabel = "Incoming interaction strength", 
                                                  title = NULL, font.size = 10, font.size.title = 10, opc_flag = F, gender_flag = F, F_flag = F, M_flag = F,
                                                  gender_flag2 = F, opc_flag2 = F, do.label = TRUE, show.legend = TRUE, show.axes = TRUE) 
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
    
    tmp <- paste0(dataset[i],".rds")
    data <- readRDS(tmp)
    
    ######################################################
    obj <- aggregateNet(data, signaling = signaling, remove.isolate = FALSE)
    num.link[[i]] <- obj@net$count
    stre.link[[i]] <- obj@net$weight
    
    num.link[[i]] <- rowSums(num.link[[i]]) + colSums(num.link[[i]]) - diag(num.link[[i]])
    stre.link[[i]] <- rowSums(stre.link[[i]]) + colSums(stre.link[[i]]) - diag(stre.link[[i]])
    
  }
  outgoing.cells <- unlist(outgoing.cells)
  incoming.cells <- unlist(incoming.cells)
  num.link <- unlist(num.link)
  stre.link <- unlist(stre.link)
  slope <- incoming.cells/outgoing.cells
  ed <- sqrt((outgoing.cells)^2+(incoming.cells)^2)
  
  group <- sub(".*--", "", names(incoming.cells))
  labels <- sub("--.*", "", names(incoming.cells))
  
  # the levels of clusters in the merged object are fixed by the levels of clusters in the female_case dataset
  df <- data.frame(x = outgoing.cells, y = incoming.cells, slope = slope, ed = ed,
                   labels = factor(labels, levels = levels(object@idents[[2]])), Interactions = num.link, Strength=stre.link, Group = factor(group, levels = unique(group)))

  write.csv(df,"signalingRole.csv")
  
  df <- df %>% mutate(Gender = factor(ifelse(str_detect(Group, "Male"), "Male", "Female")))
  df <- df %>% mutate(Condition = factor(ifelse(str_detect(Group, "MDD"), "MDD", "Control"),levels = c("MDD", "Control")))
  df <- df %>% mutate(Female = ifelse(str_detect(Group, "Female"), ifelse(str_detect(Group, "MDD"), "MDD", "Control"), "Others"))
  df <- df %>% mutate(Male = ifelse(str_detect(Group, "Male"), ifelse(str_detect(Group, "MDD"), "MDD", "Control"), "Others"))
  #df <- df %>% mutate(Group2 = ifelse(str_detect(Group, "Male"), ifelse(str_detect(Group, "MDD"), "MDD", "Control"), "Others"))
  #df$legend_group <- interaction(df$Gender, df$Condition,sep = "_")
  #df <- df %>% mutate(legend_group = paste(Gender, Condition, sep = "_")) 
  
  opc_colors <- c(Oligodendroglia = "#EEAD0E", Others = "#CCCCCC") 
  df <- df %>% mutate(Cells = ifelse(str_detect(labels, "Oli|OPC"), "Oligodendroglia", "Others"))
  opc_colors2 <- c(OPC = "#EEAD0E", OL="#DDCC77", Ast = "#A3C4DC", End = "#aaabbb", ExN="#EBA0D4", InN="#abddbc",Mic="#8C96C6",Mix="#ddcccc")
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
  
  #write.csv(df,"test.csv")
  
  gender_colors <- c(Female = "#D62728",Male = "#4682B4")
  gender_colors2 <- c(Female_MDD = "#D62728", Female_control = "#F7B7B0", Male_MDD = "#4682B4", Male_control = "#A3C4DC")
  F_case_colors <- c(MDD = "#D62728", Control = "#F7B7B0", Male = "#CCCCCC")
  M_case_colors <- c(MDD = "#4682B4", Control = "#A3C4DC", Female = "#CCCCCC") 
  
  legend_labels <- c(
    "Female_MDD" = "MDD",
    "Female_control" = "Control",
    "Male_MDD" = "MDD",
    "Male_control" = "Control"
  )
  
  gg <- ggplot(data = df, aes(x, y))
  if (opc_flag){
    gg <- gg + geom_point(aes(size = Strength, colour = Cells, fill = Cells, shape = Condition), alpha = dot.alpha) + 
      scale_fill_manual(values = opc_colors, drop = FALSE) + 
      scale_colour_manual(values = opc_colors, drop = FALSE)
  }
  else if (opc_flag2){
    gg <- gg + geom_point(aes(size = Strength, colour = Cells_Types, fill = Cells_Types, shape = Condition), alpha = dot.alpha) + 
      scale_fill_manual(values = opc_colors2, drop = FALSE) + 
      scale_colour_manual(values = opc_colors2, drop = FALSE)
  }
  else if (gender_flag){  # two colors
    gg <- gg + geom_point(aes(size = Strength, colour = Gender, fill = Gender, shape = Condition), alpha = dot.alpha) + 
      scale_fill_manual(values = gender_colors, drop = FALSE) + 
      scale_colour_manual(values = gender_colors, drop = FALSE)
  }
  else if (gender_flag2){   # four colors
    gg <- gg + geom_point(aes(size = Strength, colour = Group, fill = Group, shape = Condition), alpha = dot.alpha) + 
      scale_fill_manual(values = gender_colors2, labels=legend_labels, drop = FALSE) + 
      scale_colour_manual(name = "Group",breaks = c("Female_MDD", "Female_control", "Male_MDD", "Male_control"),values = gender_colors2, labels=legend_labels, drop = FALSE)
  }
  else if (F_flag){
    df$point_alpha <- ifelse(df$Female == "Others", dot.alpha *0.3, dot.alpha)
    gg <- gg + geom_point(aes(size = Strength, colour = Female, fill = Female, shape = Condition), alpha = df$point_alpha) + 
      scale_fill_manual(values = F_case_colors, drop = FALSE) + 
      scale_colour_manual(values = F_case_colors, drop = FALSE) + guides(colour = guide_legend(order = 1,override.aes = list(size = 3, shape = c(16, 15))),shape = "none")
  }
  else if (M_flag){
    df$point_alpha <- ifelse(df$Male == "Others", dot.alpha *0.3, dot.alpha)
    gg <- gg + geom_point(aes(size = Strength, colour = Male, fill = Male, shape = Condition), alpha = df$point_alpha) + 
      scale_fill_manual(values = M_case_colors, drop = FALSE) + 
      scale_colour_manual(values = M_case_colors, drop = FALSE)+  guides(colour = guide_legend(order = 1,override.aes = list(size = 3, shape = c(16, 15))),shape = "none")
  }
  else{
    if (is.null(color.use)) {
      color.use <- scPalette(length(unique(df$labels)))
    }
    gg <- gg + geom_point(aes(size = Strength, colour = labels, fill = labels, shape = Group), alpha = dot.alpha) + 
      scale_fill_manual(values = color.use, drop = FALSE) + 
      scale_colour_manual(values = color.use, drop = FALSE) +  guides(fill = "none", colour = "none") 
  }
  
  gg <- gg + CellChat_theme_opts() + 
    theme(text = element_text(size = font.size), 
          legend.key.height = grid::unit(0.15, "in")) + 
    labs(title = title, x = xlabel, y = ylabel) + 
    theme(plot.title = element_text(size = font.size.title, face = "plain"), 
          axis.line.x = element_line(size = 0.25), 
          axis.line.y = element_line(size = 0.25))
  
  #gg <- gg + scale_alpha_continuous(name="Strength", range = c(0.1,0.9)) + 
  #  guides(alpha = guide_legend(override.aes = list(size = 3))) 
  
  # with shapes
  if (!is.null(group)) {
    gg <- gg + scale_shape_manual(values = point.shape[1:length(unique(df$Condition))]) +  
      theme(legend.text = element_text(size = 10))
  }
  
  if (is.null(weight.MinMax)) {
    gg <- gg + scale_size_continuous(range = dot.size)
    
  } else {
    gg <- gg + scale_size_continuous(limits = weight.MinMax, range = dot.size)
  }
  
  if (do.label) {
    if (opc_flag) {
      gg <- gg + 
        ggrepel::geom_text_repel(
          aes(label = labels, colour = Cells),  
          size = label.size, show.legend = FALSE, 
          segment.size = 0.2, segment.alpha = 0.5, max.overlaps = 15
        )
      gg <- gg + guides(colour = guide_legend(order=1, override.aes = list(size = 3)),fill = "none",
                        shape = guide_legend(order=2, override.aes = list(size = 3, shape = legend.shape[1:length(unique(df$Condition))])),
                        size = guide_legend(order=3))
      gg <- gg + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, size=1, alpha = 0.05, aes(color = Cells)) + 
        scale_color_manual(values = opc_colors) +
        guides(linetype = "none") 
      
      #results <- df %>% group_by(Cells) %>% summarize(tidy(lm(y ~ x, data = cur_data()))) %>% filter(term == "x") %>% select(Cells, estimate, std.error, p.value)
    }
    else if (opc_flag2) {
      gg <- gg + 
        ggrepel::geom_text_repel(
          aes(label = labels, colour = Cells_Types),  
          size = label.size, show.legend = FALSE, 
          segment.size = 0.2, segment.alpha = 0.5, max.overlaps = 15
        )
      gg <- gg + guides(colour = guide_legend(order=1, override.aes = list(size = 3)),fill = "none",
                        shape = guide_legend(order=2, override.aes = list(size = 3, shape = legend.shape[1:length(unique(df$Condition))])),
                        size = guide_legend(order=3))
      gg <- gg + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, size=1, alpha = 0.05, aes(color = Cells_Types)) + 
        scale_color_manual(values = opc_colors2) +
        guides(linetype = "none") 
      
      #results <- df %>% group_by(Cells_Types) %>% summarize(tidy(lm(y ~ x, data = cur_data()))) %>% filter(term == "x") %>% select(Cells_Types, estimate, std.error, p.value)
    }
    else if (gender_flag) {
      gg <- gg + 
        ggrepel::geom_text_repel(
          aes(label = labels, colour = Gender),  
          size = label.size, show.legend = FALSE, 
          segment.size = 0.2, segment.alpha = 0.5, max.overlaps = 15
        )
      gg <- gg + guides(colour = guide_legend(order=1, override.aes = list(size = 3)),fill = "none",
                        shape = guide_legend(order=2, override.aes = list(size = 3, shape = legend.shape[1:length(unique(df$Condition))])),
                        size = guide_legend(order=3))
    } else if (gender_flag2) {
      gg <- gg + 
        ggrepel::geom_text_repel(
          aes(label = labels, colour = Group),  
          size = label.size, show.legend = FALSE, 
          segment.size = 0.2, segment.alpha = 0.5, max.overlaps = 15
        )
      gg <- gg + guides(colour = guide_legend(order=1, override.aes = list(shape=c(16,15,16,15),size = 3)),fill = "none", shape = F,
                        #shape = guide_legend(order=2, override.aes = list(size = 3, shape = legend.shape[1:length(unique(df$Condition))])),
                        size = guide_legend(order=2))
    }else if (F_flag) {
      gg <- gg + 
        ggrepel::geom_text_repel(
          aes(label = labels, colour = Female),  
          size = label.size, show.legend = FALSE, 
          segment.size = 0.2, segment.alpha = 0.5, max.overlaps = 15
        )
      # with shapes
      # gg <- gg + guides(colour = guide_legend(order=1, override.aes = list(size = 3)),fill = "none", 
      #                   shape = guide_legend(order=2, override.aes = list(size = 3, shape = legend.shape[1:length(unique(df$Group))])),
      #                   size = guide_legend(order=3))
      gg <- gg + guides(colour = guide_legend(order=1, override.aes = list(size = 3,shape = c(16, 15))),fill = "none",
                        # shape = guide_legend(order=2, override.aes = list(size = 3, shape = legend.shape[1:length(unique(df$Condition))])),
                        size = guide_legend(order=2))
    }
    else if (M_flag) {
      gg <- gg + 
        ggrepel::geom_text_repel(
          aes(label = labels, colour = Male),  
          size = label.size, show.legend = FALSE, 
          segment.size = 0.2, segment.alpha = 0.5, max.overlaps = 15
        )
      gg <- gg + guides(colour = guide_legend(order=1, override.aes = list(size = 3,shape = c(16, 15))),fill = "none",
                        # shape = guide_legend(order=2, override.aes = list(size = 3, shape = legend.shape[1:length(unique(df$Condition))])),
                        size = guide_legend(order=2))
    }
    else {
      
      gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels, colour = labels), 
                                          size = label.size, show.legend = FALSE, 
                                          segment.size = 0.2, segment.alpha = 0.5,max.overlaps = 15)
      gg <- gg + guides(shape = guide_legend(order=1, override.aes = list(size = 3, shape = legend.shape[1:length(unique(df$Group))])),
                        size = guide_legend(order=2)) # no color legend
    }
  }
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }
  
  if (!show.axes) {
    gg <- gg + theme_void()
  }
  print(gg)
}
