netAnalysis_signalingChanges_scatter_V2 <- function (object, idents.use.multi, color.use = c("grey10", "#F8766D", 
                                            "#00BFC4"), comparison = c(1, 2), signaling = NULL, signaling.label = NULL, 
          top.label = 1, signaling.exclude = NULL, xlims = NULL, ylims = NULL, 
          slot.name = "netP", dot.size = 2.5, point.shape = c(16,15,17,18), label.size = 3, dot.alpha = 0.7, x.measure = "outdeg", 
          y.measure = "indeg", xlabel = "Differential outgoing interaction strength", 
          ylabel = "Differential incoming interaction strength", title = NULL, 
          font.size = 10, font.size.title = 10, do.label = T, show.legend = T, 
          show.axes = T) 
{
  #idents.use <- unique(cellchat@idents[[2]])[41]
  #print(idents.use.multi)
  #======================================multiple cluster section==========================
  if (length(strsplit(idents.use.multi, split=" ")) >=1) {
    message(paste0("Visualizing differential outgoing and incoming signaling changes of ", length(strsplit(idents.use.multi, split=" "))," Clusters..."))
    idents.use.all <- unlist(strsplit(idents.use.multi, split = " "))
    
    if (is.list(object)) {
      object <- mergeCellChat(object, add.names = names(object))
    }
    if (is.list(object@net[[1]])) {
      dataset.name <- names(object@net)
      message(paste0("Visualizing differential outgoing and incoming signaling changes from ", 
                     dataset.name[comparison[1]], " to ", dataset.name[comparison[2]]))
      #title <- paste0("Signaling changes of ", idents.use.all[1], idents.use.all[2],idents.use.all[3],
      #                " (", dataset.name[comparison[1]], " vs. ", dataset.name[comparison[2]], 
      #                ")")
      title <- paste0("Signaling changes of ", dataset.name[comparison[1]], " vs. ", dataset.name[comparison[2]])
      
      cell.levels <- levels(object@idents$joint)
      if (is.null(xlabel) | is.null(ylabel)) {
        xlabel = "Differential outgoing interaction strength"
        ylabel = "Differential incoming interaction strength"
      }
    }
    else {
      message("Visualizing outgoing and incoming signaling on a single object \n")
      title <- paste0("Signaling patterns of ", idents.use.all)
      if (length(slot(object, slot.name)$centr) == 0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
      }
      cell.levels <- levels(object@idents)
    }
    
    # ************************************** parsing individual cluster *************************************
    df <- data.frame(outgoing=numeric(),incoming=numeric(),specificity.out.in=numeric(),specificity=numeric(),cluster=character(),stringsAsFactors = F)

    out.ratio.all <- data.frame(outgoing=numeric(),stringsAsFactors = F)
    in.ratio.all <- data.frame(incoming=numeric(),stringsAsFactors = F)
    for (ij in 1:length(idents.use.all)){

    idents.use <- idents.use.all[ij]
    #print(idents.use)

    if (!(idents.use %in% cell.levels)) {
      stop("Please check the input cell group names!")
    }
    if (is.null(signaling)) {
      signaling <- union(object@netP[[comparison[1]]]$pathways, 
                         object@netP[[comparison[2]]]$pathways)
    }
    if (!is.null(signaling.exclude)) {
      signaling <- setdiff(signaling, signaling.exclude)
    }
    mat.all.merged <- list()
    for (ii in 1:length(comparison)) {
      if (length(slot(object, slot.name)[[comparison[ii]]]$centr) == 
          0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores for each dataset seperately! ")
      }
      
      # check whether the "outdeg" and "indeg" columns are represented in the names of centr list of one certain dataset
      if (sum(c(x.measure, y.measure) %in% names(slot(object, 
                                                      slot.name)[[comparison[ii]]]$centr[[1]])) != 2) {
        stop(paste0("`x.measure, y.measure` should be one of ", 
                    paste(names(slot(object, slot.name)[[comparison[ii]]]$centr[[1]]), 
                          collapse = ", "), "\n", "`outdeg_unweighted` is only supported for version >= 1.1.2"))
      }
      centr <- slot(object, slot.name)[[comparison[ii]]]$centr
      outgoing <- matrix(0, nrow = length(cell.levels), ncol = length(centr))
      incoming <- matrix(0, nrow = length(cell.levels), ncol = length(centr))
      dimnames(outgoing) <- list(cell.levels, names(centr))
      dimnames(incoming) <- dimnames(outgoing)
      for (i in 1:length(centr)) {
        outgoing[, i] <- centr[[i]][[x.measure]]
        incoming[, i] <- centr[[i]][[y.measure]]
      }
      mat.out <- t(outgoing)  # mat.out: pathways of one certain dataset*cluster
      mat.in <- t(incoming)
      mat.all <- array(0, dim = c(length(signaling), ncol(mat.out), 
                                  2))
      mat.t <- list(mat.out, mat.in)
      
      #parsing the outgoing and incoming matrix separately
      for (i in 1:length(comparison)) {
        mat = mat.t[[i]]      # mat: pathway*cluster 
        mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]      # mat1: matched pathways*all clusters, ensuring the format of matrix by DROP=FALSE. matched: overlapped pathways between this dataset and all "union/setdiff" pathways
        mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))   # new mat: a larger matrix of union/setdiff pathways*cluster
        idx <- match(rownames(mat1), signaling)                        # idx: row index for new mat, length(idx) <= No. pathways of one certain dataset
        mat[idx[!is.na(idx)], ] <- mat1                                # assignment of the larger matrix with matched pathways of one certain dataset
        dimnames(mat) <- list(signaling, colnames(mat1))
        mat.all[, , i] = mat
      }
      # mat.all:[[signaling,cluster,outgoing],[signaling,cluster,incoming]]
      dimnames(mat.all) <- list(dimnames(mat)[[1]], dimnames(mat)[[2]], c("outgoing", "incoming"))
      
      # a list of two datasets with parsed outgoing and incoming matrix
      # mat.all.merged: [[1]] dataset1, [[1]] outgoing, [[2]] incoming; [[2]] dataset2, [[1]] outgoing, [[2]] incoming.
      mat.all.merged[[ii]] <- mat.all
    }
    mat.all.merged.use <- list(mat.all.merged[[1]][, idents.use, ], mat.all.merged[[2]][, idents.use, ])
    idx.specific <- mat.all.merged.use[[1]] * mat.all.merged.use[[2]]   # multiplication of corresponding elements in the two tensors
    #print("..................")
    #print(idx.specific)           #      for finding the '0' signals in two datasets (either '0' works)
    mat.sum <- mat.all.merged.use[[2]] + mat.all.merged.use[[1]]        # sum of corresponding elements in the two tensors
    #print(",,,,,,,,,,,,,,,,,,")
    #print(mat.sum)                #      for finding the 'both 0' signals in two datsets
    out.specific.signaling <- rownames(idx.specific)[(mat.sum[, 1] != 0) & (idx.specific[, 1] == 0)]   # outgoing signals existing 
    in.specific.signaling <- rownames(idx.specific)[(mat.sum[, 2] != 0) & (idx.specific[, 2] == 0)]    # incoming signals existing
    nonzero.out.signaling <- rownames(idx.specific)[idx.specific[, 1] != 0]
    nonzero.in.signaling <- rownames(idx.specific)[idx.specific[, 2] != 0]
    #dif.fold <- mat.all.merged.use[[1]] / mat.all.merged.use[[2]]
    nonzero.out.mat1 <- mat.all.merged.use[[1]][rownames(mat.all.merged.use[[1]]) %in% nonzero.out.signaling,]
    nonzero.out.mat2 <- mat.all.merged.use[[2]][rownames(mat.all.merged.use[[2]]) %in% nonzero.out.signaling,]
    nonzero.in.mat1 <- mat.all.merged.use[[1]][rownames(mat.all.merged.use[[1]]) %in% nonzero.in.signaling,]
    nonzero.in.mat2 <- mat.all.merged.use[[2]][rownames(mat.all.merged.use[[2]]) %in% nonzero.in.signaling,]
    
    #selected_mat <- mat[rownames(mat) %in% rows_to_keep, ]
    #print(class(nonzero.out.mat1/nonzero.out.mat2))
    out.ratio <- nonzero.out.mat1/nonzero.out.mat2
    in.ratio <- nonzero.in.mat1/nonzero.in.mat2
    if (is.null(dim(in.ratio))){
      in.ratio <- t(as.matrix(in.ratio))
      rownames(in.ratio) <-  nonzero.in.signaling
      in.ratio.des <- in.ratio
    }
    else
    {
      in.ratio.des <- in.ratio[order(in.ratio[,2],decreasing = T),]
    }
    if (is.null(dim(out.ratio))){
      out.ratio <- t(as.matrix(out.ratio))
      rownames(out.ratio) <-  nonzero.out.signaling
      out.ratio.des <- out.ratio
    }
    else
    {
    out.ratio.des <- out.ratio[order(out.ratio[,1],decreasing = T),]
    }
    #print(class(nonzero.in.mat1/nonzero.in.mat2))
    #print(nonzero.in.mat1/nonzero.in.mat2)
    #colnames(out.ratio.des)[1] <- paste0(idents.use,"_outgoing")
    out.ratio.des <- as.data.frame(out.ratio.des)
    out.ratio.des$Cluster <- idents.use
    out.ratio.all <- rbind(out.ratio.all, out.ratio.des)

    in.ratio.des <- as.data.frame(in.ratio.des)
    in.ratio.des$Cluster <- idents.use
    in.ratio.all <- rbind(in.ratio.all, in.ratio.des)

    mat.diff <- mat.all.merged.use[[2]] - mat.all.merged.use[[1]]
    
    #print(mat.diff)
    # all the idx.specific, mat.sum, mat.diff look like as follows, as the xxx comes from the operations between two datasets
    #            outgoing     incoming
    # signals      xxx           xxx
    
    idx <- rowSums(mat.diff) != 0               # summing up the signals of outgoing([,1]) and incoming([,2])
    #print(idx)
    mat.diff <- mat.diff[idx, ]
    #print(dim(mat.diff))
    out.specific.signaling <- rownames(mat.diff) %in% out.specific.signaling   # outgoing signals existing differently in two datasets
    in.specific.signaling <- rownames(mat.diff) %in% in.specific.signaling     # incoming signals existing differently in two datasets
    
    out.in.specific.signaling <- as.logical(out.specific.signaling *     # outgoing and incoming signals both existing differently in two datasets
                                              in.specific.signaling)
    specificity.out.in <- matrix(0, nrow = nrow(mat.diff), ncol = 1)
    specificity.out.in[out.in.specific.signaling] <- 2                   # shared signals <- 2
    specificity.out.in[setdiff(which(out.specific.signaling), which(out.in.specific.signaling))] <- 1    # outgoing specific signals <- 1
    specificity.out.in[setdiff(which(in.specific.signaling), which(out.in.specific.signaling))] <- -1    # incoming specific signals <- -1
    
    specificity <- numeric(dim(specificity.out.in)[1])
    specificity[(specificity.out.in != 0) & (rowSums(mat.diff >= 0) == 2)] = 1
    specificity[(specificity.out.in != 0) & (rowSums(mat.diff <= 0) == 2)] = -1
    #df <- as.data.frame(mat.diff)
    #df$specificity.out.in <- specificity.out.in
    #df$specificity = 0
    #df$specificity[(specificity.out.in != 0) & (rowSums(mat.diff >= 0) == 2)] = 1
    #df$specificity[(specificity.out.in != 0) & (rowSums(mat.diff <= 0) == 2)] = -1
    #df$cluster = idents.use

    df <- rbind(df, data.frame(outgoing=mat.diff[,1],incoming=mat.diff[,2],specificity.out.in=specificity.out.in,specificity=specificity,cluster=idents.use))

    }
    #write.csv(out.ratio.all,"out.ratio.all.csv")
    #write.csv(in.ratio.all,"in.ratio.all.csv")

  #==================================All clusters parsed==============================================
  
  out.in.category <- c("Shared", "Incoming specific", "Outgoing specific", "Incoming & Outgoing specific")
  specificity.category <- c("Shared", paste0(dataset.name[comparison[1]], " specific"), paste0(dataset.name[comparison[2]], " specific"))
  cluster.category <- c("Ast1", "Ast2", "End1", "ExN1_L24","ExN10_L46","ExN11_L56","ExN12_L56","ExN13_L56","ExN14","ExN15_L56","ExN16_L56", "ExN17","ExN18","ExN19_L56","ExN2_L23","ExN20_L56","ExN3_L46","ExN4_L35","ExN5","ExN6","ExN7", "ExN8_L24","ExN9_L23","InN1_PV","InN10_ADARB2","InN2_SST","InN3_VIP","InN4_VIP","InN5_SST","InN6_LAMP5","InN7_Mix","InN8_Mix","InN9_PV","Mic1","Mix","Oli1","Oli2","Oli3","OPC1","OPC2","OPC3")
  cluster.category <- factor(cluster.category, levels = cluster.category)
  col_cluster.category <- scPalette(length(cluster.category))
  #col_cluster.category <- factor(col_cluster.category,levels = col_cluster.category)
  color_mapping <- setNames(col_cluster.category,cluster.category)

  df$specificity.out.in <- plyr::mapvalues(df$specificity.out.in, from = c(0, -1, 1, 2), to = out.in.category)
  df$specificity.out.in <- factor(df$specificity.out.in, levels = out.in.category)
  df$specificity <- plyr::mapvalues(df$specificity, from = c(0, -1, 1), to = specificity.category)
  df$specificity <- factor(df$specificity, levels = specificity.category)
  df$cluster <- factor(df$cluster,levels=cluster.category)

  # shapes correspond to datasets
  #point.shape.use <- point.shape[out.in.category %in% unique(df$specificity.out.in)]
  point.shape.use <- point.shape[specificity.category %in% unique(df$specificity)]

  df$specificity.out.in = droplevels(df$specificity.out.in, exclude = setdiff(out.in.category, unique(df$specificity.out.in)))
  df$specificity = droplevels(df$specificity, exclude = setdiff(specificity.category, unique(df$specificity)))

  # colors correspond to clusters
  df$cluster <- droplevels(df$cluster, exclude = setdiff(cluster.category, unique(df$cluster)))

  #color.use <- color.use[specificity.category %in% unique(df$specificity)]
  #color.use <- col_cluster.category[match(idents.use.all, cluster.category)]
  #color.use <- droplevels(color.use, exclude = setdiff(col_cluster.category, color.use))

  color.use <- color_mapping[idents.use.all]
  df$labels <- rownames(df)
  }


  #==========================================mapping====================================
  
  #gg <- ggplot(data = df, aes(outgoing, incoming)) + geom_point(aes(colour = specificity, 
  #                                                                  fill = specificity, shape = specificity.out.in), size = dot.size)
  # write.csv(df,"F2_local_ChangeScatter.csv")
  gg <- ggplot(data = df, aes(outgoing, incoming)) + geom_point(aes(colour = cluster, 
                                                                  fill = cluster, shape = specificity), size = dot.size)
  gg <- gg + theme_linedraw() + theme(panel.grid = element_blank()) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", 
               size = 0.25) + geom_vline(xintercept = 0, linetype = "dashed", 
                                         color = "grey50", size = 0.25) + theme(text = element_text(size = font.size), 
                                                                                legend.key.height = grid::unit(0.15, "in")) + labs(title = title, 
                                                                                                                                   x = xlabel, y = ylabel) + theme(plot.title = element_text(size = font.size.title, 
                                                                                                                                                                                             hjust = 0.5, face = "plain")) + theme(axis.line.x = element_line(size = 0.25), 
                                                                                                                                                                                                                                   axis.line.y = element_line(size = 0.25))
  gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use, 
                                                       alpha = dot.alpha), drop = FALSE) + guides(fill = "none")
  gg <- gg + scale_colour_manual(values = color.use, drop = FALSE)
  gg <- gg + scale_shape_manual(values = point.shape.use)
  gg <- gg + theme(legend.title = element_blank())
  if (!is.null(xlims)) {
    gg <- gg + xlim(xlims)
  }
  if (!is.null(ylims)) {
    gg <- gg + ylim(ylims)
  }
  if (do.label) {
    if (is.null(signaling.label)) {
      thresh <- stats::quantile(abs(as.matrix(df[, 1:2])), 
                                probs = 1 - top.label)
      idx = abs(df[, 1]) > thresh | abs(df[, 2]) > thresh
      data.label <- df[idx, ]
    }
    else {
      data.label <- df[rownames(df) %in% signaling.label, 
      ]
    }
    gg <- gg + ggrepel::geom_text_repel(data = data.label, 
                                        mapping = aes(label = labels, colour = cluster), 
                                        size = label.size, show.legend = F, segment.size = 0.2, 
                                        segment.alpha = 0.5)
  }
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }
  if (!show.axes) {
    gg <- gg + theme_void()
  }
  gg
}