netAnalysis_contribution_v2 <- function (object, dataset = c(NULL,NULL), signaling, signaling.name = NULL, sources.use = NULL, 
    targets.use = NULL, width = 0.1, vertex.receiver = NULL, 
    thresh = 0.05, return.data = FALSE, x.rotation = 0, title = "Contribution of each L-R pair", 
    font.size = 10, font.size.title = 10) 
{
  if (length(dataset) == 2) {    # two datasets
    df.use <- data.frame(name=character(),contribution=numeric(),data=character(),stringsAsFactors = F)
    pSum.use <- list()
    for (i in 1:length(dataset)) {
      pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR[[i]]$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)
      pair.name.use = dplyr::select(object@DB$interaction[rownames(pairLR), ], "interaction_name_2")
      if (is.null(signaling.name)) {
        signaling.name <- signaling
      }
      net <- object@net[[dataset[i]]]
      pairLR.use.name <- dimnames(net$prob)[[3]]
      pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
      pairLR <- pairLR[pairLR.name, ]
      prob <- net$prob
      pval <- net$pval
      prob[pval > thresh] <- 0
      if (!is.null(sources.use)) {
        if (is.character(sources.use)) {
          if (all(sources.use %in% dimnames(prob)[[1]])) {
            sources.use <- match(sources.use, dimnames(prob)[[1]])
          }
          else {
            stop("The input `sources.use` should be cell group names or a numerical vector!")
          }
        }
        idx.t <- setdiff(1:nrow(prob), sources.use)
        prob[idx.t, , ] <- 0
      }
      if (!is.null(targets.use)) {
        if (is.character(targets.use)) {
          if (all(targets.use %in% dimnames(prob)[[1]])) {
            targets.use <- match(targets.use, dimnames(prob)[[2]])
          }
          else {
            stop("The input `targets.use` should be cell group names or a numerical vector!")
          }
        }
        idx.t <- setdiff(1:nrow(prob), targets.use)
        prob[, idx.t, ] <- 0
      }
      if (length(pairLR.name) > 1) {
        pairLR.name.use <- pairLR.name[apply(prob[, , pairLR.name], 3, sum) != 0]
      }
      else {
        pairLR.name.use <- pairLR.name[sum(prob[, , pairLR.name]) != 0]
      }
      if (length(pairLR.name.use) == 0) {
        stop(paste0("There is no significant communication of ", 
                    signaling.name))
      }
      else {
        pairLR <- pairLR[pairLR.name.use, ]
      }
      prob <- prob[, , pairLR.name.use]  # "There are significant communications."
      if (length(dim(prob)) == 2) {
        prob <- replicate(1, prob, simplify = "array")
        dimnames(prob)[3] <- pairLR.name.use
      }
      prob1 <- prob 
      prob <- (prob - min(prob))/(max(prob) - min(prob))   # Min-Max Normalization
#      prob.use[[i]] <- prob
      
      if (is.null(vertex.receiver)) {
        pSum <- apply(prob, 3, sum)
        pSum1 <- apply(prob1, 3, sum)
        pSum.max <- sum(prob)
        pSum.max1 <- sum(prob1)

        pSum <- pSum/pSum.max
        pSum[is.na(pSum)] <- 0
        #pSum1 <- pSum1/pSum.max1
        #pSum1[is.na(pSum1)] <- 0
        #y.lim <- max(pSum)
        pair.name <- unlist(dimnames(prob)[3])
        pair.name <- factor(pair.name, levels = unique(pair.name))
        if (!is.null(pairLR.name.use)) {
          pair.name <- pair.name.use[as.character(pair.name), 1]
          pair.name <- factor(pair.name, levels = unique(pair.name))
        }
      }
      pSum.use[[i]] <- pSum1    # non-scaled values
    }
        mat <- pSum.use[[2]] - pSum.use[[1]]   # percentage
        y.lim <- max(abs(mat))
        df1 <- data.frame(name = pair.name, contribution = abs(mat))
        
          df <- df1
        
        df <- df[order(df$contribution, decreasing = TRUE), ]
        df$name <- factor(df$name, levels = df$name[order(df$contribution, decreasing = TRUE)])
        df1$name <- factor(df1$name, levels = df1$name[order(df1$contribution, decreasing = TRUE)])
        print(df)
        
        nrxn3_nlgn1_rows <- grep("^NRXN3 - NLGN1", df$name, value = TRUE, ignore.case = TRUE)
        nrxn1_nlgn1_rows <- grep("^NRXN1 - NLGN1", df$name, value = TRUE, ignore.case = TRUE)
        nrxn3_nlgn1_data <- df[df$name %in% nrxn3_nlgn1_rows, , drop = FALSE]
        nrxn1_nlgn1_data <- df[df$name %in% nrxn1_nlgn1_rows, , drop = FALSE]
        other_data <- df[!df$name %in% c(nrxn3_nlgn1_rows,nrxn1_nlgn1_rows), , drop = FALSE]

        wb <- createWorkbook()
        
        addWorksheet(wb, "Pairs")
        addWorksheet(wb, "NRXN3_NLGN1")
        addWorksheet(wb, "NRXN1_NLGN1")
        addWorksheet(wb, "Others")
        
        writeData(wb, sheet = "Pairs", df)
        writeData(wb, sheet = "NRXN3_NLGN1", nrxn3_nlgn1_data)
        writeData(wb, sheet = "NRXN1_NLGN1", nrxn1_nlgn1_data)
        writeData(wb, sheet = "Others", other_data)
        
        # saveWorkbook(wb, "D:/Project/niu/OPC-MDD/permutation/Pairs.xlsx", overwrite = TRUE)
        saveWorkbook(wb, "Pairs.xlsx", overwrite = TRUE)
        
        
  }
}
