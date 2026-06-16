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
      mat.all <- array(0, dim = c(length(signaling), ncol(mat.out), 2))
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

    df <- rbind(df, data.frame(outgoing=mat.diff[,1],incoming=mat.diff[,2],specificity.out.in=specificity.out.in,specificity=specificity,cluster=idents.use))
    }
    library(tibble)
    df <- rownames_to_column(df, var = "pathway")
    # write.csv(df,"D:/Project/niu/OPC-MDD/permutation/df.csv", row.names = FALSE)
    write.csv(df,"df.csv", row.names = FALSE)
    
    df_pathway <- df %>% select(pathway,outgoing,incoming,cluster) 
    colnames(df_pathway) <- c("pathway","x","y","cell_type")
    
    library(dplyr)
    library(tidyr)
    library(openxlsx)
    
    df_pathway <- df_pathway %>%
      mutate(
        ed = sqrt((x)^2+(y)^2)
      )
    # write.csv(df_pathway,"D:/Project/niu/OPC-MDD/permutation/df_pathway.csv", row.names = FALSE)
    write.csv(df_pathway,"df_pathway.csv", row.names = FALSE)
    
    nrxn_rows <- grep("^NRXN", df_pathway$pathway, value = TRUE, ignore.case = TRUE)
    nrxn_data <- df_pathway[df_pathway$pathway %in% nrxn_rows, , drop = FALSE]
    
    other_data <- df_pathway[!df_pathway$pathway %in% nrxn_rows, , drop = FALSE]
    
    wb <- createWorkbook()
    
    addWorksheet(wb, "AllPathways")
    addWorksheet(wb, "NRXN")
    addWorksheet(wb, "Others")
    
    writeData(wb, sheet = "AllPathways", df_pathway)
    writeData(wb, sheet = "NRXN", nrxn_data)
    writeData(wb, sheet = "Others", other_data)
    
    # saveWorkbook(wb, "D:/Project/niu/OPC-MDD/permutation/Pathways.xlsx", overwrite = TRUE)
    saveWorkbook(wb, "Pathways.xlsx", overwrite = TRUE)
    
  }
}