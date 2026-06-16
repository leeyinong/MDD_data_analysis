netAnalysis_contribution_v2 <- function (object, dataset = c(NULL,NULL), signaling, signaling.name = NULL, sources.use = NULL, 
    targets.use = NULL, width = 0.1, vertex.receiver = NULL, 
    thresh = 0.05, return.data = FALSE, x.rotation = 0, title = "Contribution of each L-R pair", 
    font.size = 10, font.size.title = 10) 
{
  print(dataset)
  if (length(dataset) == 1){
    pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)
    pair.name.use = dplyr::select(object@DB$interaction[rownames(pairLR), ], "interaction_name_2")
    if (is.null(signaling.name)) {
        signaling.name <- signaling
    }
    net <- object@net
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
    prob <- (prob - min(prob))/(max(prob) - min(prob))
    if (is.null(vertex.receiver)) {
        pSum <- apply(prob, 3, sum)
        pSum.max <- sum(prob)
        pSum <- pSum/pSum.max
        pSum[is.na(pSum)] <- 0
        y.lim <- max(pSum)
        pair.name <- unlist(dimnames(prob)[3])
        pair.name <- factor(pair.name, levels = unique(pair.name))
        if (!is.null(pairLR.name.use)) {
            pair.name <- pair.name.use[as.character(pair.name), 1]
            pair.name <- factor(pair.name, levels = unique(pair.name))
        }
        mat <- pSum   # percentage
        df1 <- data.frame(name = pair.name, contribution = mat)
        if (nrow(df1) < 10) {
            df2 <- data.frame(name = as.character(1:(10 - nrow(df1))), 
                contribution = rep(0, 10 - nrow(df1)))
            df <- rbind(df1, df2)
        }
        else {
            df <- df1
        }
        df <- df[order(df$contribution, decreasing = TRUE), ]
        df$name <- factor(df$name, levels = df$name[order(df$contribution, decreasing = TRUE)])
        df1$name <- factor(df1$name, levels = df1$name[order(df1$contribution, decreasing = TRUE)])
        print(df)
        gg <- ggplot(df, aes(x = name, y = contribution)) + geom_bar(stat = "identity", 
            width = 0.7) + theme_classic() + theme(axis.text.y = element_text(angle = x.rotation, 
            hjust = 1, size = font.size, colour = "black"), axis.text = element_text(size = font.size), 
            axis.title.y = element_text(size = font.size), axis.text.x = element_blank(), 
            axis.ticks = element_blank()) + xlab("") + ylab("Relative contribution") + 
            ylim(0, y.lim) + coord_flip() + theme(legend.position = "none") + 
            scale_x_discrete(limits = rev(levels(df$name)), labels = c(rep("", 
                max(0, 10 - nlevels(df1$name))), rev(levels(df1$name))))
        if (!is.null(title)) {
            gg <- gg + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, 
                size = font.size.title))
        }
        gg
    }
    else {
        pair.name <- factor(unlist(dimnames(prob)[3]), levels = unique(unlist(dimnames(prob)[3])))
        pSum <- apply(prob, 3, sum)
        pSum.max <- sum(prob)
        pSum <- pSum/pSum.max
        pSum[is.na(pSum)] <- 0
        y.lim <- max(pSum)
        df <- data.frame(name = pair.name, contribution = pSum)
        gg <- ggplot(df, aes(x = name, y = contribution)) + geom_bar(stat = "identity", 
            width = 0.2) + theme_classic() + theme(axis.text = element_text(size = 10), 
            axis.text.x = element_text(angle = x.rotation, hjust = 1, 
                size = 8), axis.title.y = element_text(size = 10)) + 
            xlab("") + ylab("Relative contribution") + ylim(0, 
            y.lim) + ggtitle("All") + theme(plot.title = element_text(hjust = 0.5))
        if (dim(prob)[3] > 1) {
            pSum <- apply(prob[, vertex.receiver, ], 3, sum)
        }
        else {
            pSum <- sum(prob[, vertex.receiver, ])
        }
        pSum <- pSum/pSum.max
        pSum[is.na(pSum)] <- 0
        df <- data.frame(name = pair.name, contribution = pSum)
        gg1 <- ggplot(df, aes(x = name, y = contribution)) + 
            geom_bar(stat = "identity", width = 0.2) + theme_classic() + 
            theme(axis.text = element_text(size = 10), axis.text.x = element_text(angle = x.rotation, 
                hjust = 1, size = 8), axis.title.y = element_text(size = 10)) + 
            xlab("") + ylab("Relative contribution") + ylim(0, 
            y.lim) + ggtitle("Hierarchy1") + theme(plot.title = element_text(hjust = 0.5))
        if (dim(prob)[3] > 1) {
            pSum <- apply(prob[, setdiff(1:dim(prob)[1], vertex.receiver), 
                ], 3, sum)
        }
        else {
            pSum <- sum(prob[, setdiff(1:dim(prob)[1], vertex.receiver), 
                ])
        }
        pSum <- pSum/pSum.max
        pSum[is.na(pSum)] <- 0
        df <- data.frame(name = pair.name, contribution = pSum)
        gg2 <- ggplot(df, aes(x = name, y = contribution)) + 
            geom_bar(stat = "identity", width = 0.9) + theme_classic() + 
            theme(axis.text = element_text(size = 10), axis.text.x = element_text(angle = x.rotation, 
                hjust = 1, size = 8), axis.title.y = element_text(size = 10)) + 
            xlab("") + ylab("Relative contribution") + ylim(0, 
            y.lim) + ggtitle("Hierarchy2") + theme(plot.title = element_text(hjust = 0.5))
        title <- cowplot::ggdraw() + cowplot::draw_label(paste0("Contribution of each signaling in ", 
            signaling.name, " pathway"), fontface = "bold", size = 10)
        gg.combined <- cowplot::plot_grid(gg, gg1, gg2, nrow = 1)
        gg.combined <- cowplot::plot_grid(title, gg.combined, 
            ncol = 1, rel_heights = c(0.1, 1))
        gg <- gg.combined
        gg
    }
    if (return.data) {
        df <- subset(df, contribution > 0)
        return(list(LR.contribution = df, gg.obj = gg))
    }
    else {
        return(gg)
    }
  }
  else if (length(dataset) == 2) {    # two datasets
    df.use <- data.frame(name=character(),contribution=numeric(),data=character(),stringsAsFactors = F)
    gender_colors <- c(Male = "steelblue", Female = "#d62728")
    pSum.use <- list()
    for (i in 1:length(dataset)) {
      print(dataset[i])
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
        print("pSum and pSum1")
        print(pSum)
        print(pSum1)
        
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
      pSum.use[[i]] <- pSum1
      print(pSum.use[[i]])
    }
        mat <- pSum.use[[2]] - pSum.use[[1]]   # percentage
        y.lim <- max(abs(mat))
        df1 <- data.frame(name = pair.name, contribution = abs(mat))
        if (nrow(df1) < 10) {
          df2 <- data.frame(name = as.character(1:(10 - nrow(df1))), 
                            contribution = rep(0, 10 - nrow(df1)))
          df <- rbind(df1, df2)
        }
        else {
          df <- df1
        }
        df <- df[order(df$contribution, decreasing = TRUE), ]
        df$name <- factor(df$name, levels = df$name[order(df$contribution, decreasing = TRUE)])
        df1$name <- factor(df1$name, levels = df1$name[order(df1$contribution, decreasing = TRUE)])
        print(df)
        gg <- ggplot(df, aes(x = name, y = contribution)) + geom_bar(stat = "identity", 
                                                                     width = 0.7) + theme_classic() + theme(axis.text.y = element_text(angle = x.rotation, 
                                                                                                                                       hjust = 1, size = font.size, colour = "black"), axis.text = element_text(size = font.size), 
                                                                                                            axis.title.y = element_text(size = font.size), axis.text.x = element_blank(), 
                                                                                                            axis.ticks = element_blank()) + xlab("") + ylab("Contribution") + scale_y_continuous(breaks = seq(0, 100, 10))+
          ylim(0, 100) + coord_flip() + theme(legend.position = "none") + 
          scale_x_discrete(limits = rev(levels(df$name)), labels = c(rep("", 
                                                                         max(0, 10 - nlevels(df1$name))), rev(levels(df1$name))))
        if (!is.null(title)) {
          gg <- gg + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, 
                                                                      size = font.size.title))
        }
        gg
        
      if (return.data) {
        df <- subset(df, contribution > 0)
        return(list(LR.contribution = df, gg.obj = gg))
      }
      else {
        return(gg)
      }
  }
  else {    # four datasets
    df.use <- data.frame(name=character(),contribution=numeric(),data=character(),stringsAsFactors = F)
    gender_colors <- c(Female = "#d62728",Male = "steelblue")
    gender <- c("Female","Male")
    pSum.use <- list()
    mat <- list()
    df_list <- list()
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
      if (length(pairLR.name) == 1) {
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
      pSum.use[[i]] <- pSum1
     # print(pSum.use[[i]])
    }
  for (j in 1:length(gender)) {
    mat[[j]] <- pSum.use[[2*j]] - pSum.use[[2*j-1]]   # absolute value
    y.lim <- max(abs(mat[[j]]))
    df1 <- data.frame(name = pair.name, contribution = abs(mat[[j]]), Gender = gender[j])
    
      df <- df1
    
    df_list[[j]] <- df
  }
    combined_df <- do.call(rbind,df_list)
    
    female_order <- combined_df %>% filter(Gender== "Female") %>% arrange(desc(contribution)) %>% pull(name) %>% unique()
    female_order <- factor(female_order,levels = female_order)
    combined_df$Gender <- factor(combined_df$Gender, levels = c("Female", "Male"))
    print(combined_df)
    #df <- df[order(df$contribution, decreasing = TRUE), ]
    df$name <- factor(df$name, levels = female_order)
    df1$name <- factor(df1$name, levels = female_order)
    gg <- ggplot(combined_df, aes(x = name, y = contribution, fill = Gender)) + geom_bar(stat = "identity",position = position_dodge2(width=0.7,reverse=T), 
                                                                 width = 0.7) +scale_fill_manual(values = gender_colors) + theme_classic() + theme(axis.text.y = element_text(angle = x.rotation, 
                                                                                                                                                                        hjust = 1, size = font.size, colour = "black"), axis.text = element_text(size = font.size), 
                                                                                                                                             axis.title.y = element_text(size = font.size), axis.text.x = element_blank(), 
                                                                                                                                             axis.ticks = element_blank()) + xlab("") + ylab("Strength") + scale_y_continuous(breaks = seq(0, 100, 10))+
      ylim(0, 90) + coord_flip() + theme(legend.position = "right") + 
      scale_x_discrete(limits = rev(levels(df$name))) + labs(title = "Signal alternations of each L-R pair") + theme(plot.title = element_text(hjust = 0.5))
    
    gg
    
    if (return.data) {
      df <- subset(df, contribution > 0)
      return(list(LR.contribution = df, gg.obj = gg))
    }
    else {
      return(gg)
    }
  }
}
