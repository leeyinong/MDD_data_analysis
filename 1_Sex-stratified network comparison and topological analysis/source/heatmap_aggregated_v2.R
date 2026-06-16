heatmap_aggregated_v2 <- function (object, dataset = c(NULL,NULL), method = c("weight", "count", "weighted_count", 
                             "weighted_count2", "weight_threshold"), cut_off = 0.05, interaction_use = "all", 
          group = NULL, sender.names = NULL, receiver.names = NULL) 
{
  method <- match.arg(method)
  if (length(dataset) == 1){
  if (is.null(sender.names)) {
    sender.names = rownames(object@net[[1]])
  }
  if (is.null(receiver.names)) {
    receiver.names = colnames(object@net[[1]])
  }
  if (interaction_use == "all") {
    net_aggregated <- net_aggregation_v2(object@net, method = method, 
                                      cut_off = cut_off)
  }
  else {
    net_aggregated <- net_aggregation_v2(object@net[interaction_use], 
                                      method = method, cut_off = cut_off)
  }
  ylgnbu_colors <- brewer.pal(n = 9, name = "YlGnBu")
  col_map = colorRamp2(c(0, max(net_aggregated[sender.names, 
                                               receiver.names])/2, max(net_aggregated[sender.names, 
                                                                                      receiver.names])), c("#2461A1", "#F5F9FC", "#A91E2C"))
  if (is.null(group)) {
    ComplexHeatmap::Heatmap(net_aggregated[sender.names, 
                                           receiver.names], name = method, cluster_rows = F, 
                            cluster_columns = F, column_names_rot = 45, row_names_side = "left", 
                            col = col_map, column_title = "Receiver", row_title = "Sender", 
                            column_title_side = "bottom", heatmap_legend_param = list(color_bar = "continuous"))
  }
  else {
    left_Annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = scPalette(length(unique(group[sender.names])))), 
                                                     labels = sort(unique(group[sender.names]))))
    bottom_Annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scPalette(length(unique(group[receiver.names])))), 
                                                           labels = sort(unique(group[receiver.names]))))
    ComplexHeatmap::Heatmap(net_aggregated[sender.names, 
                                           receiver.names], name = "Commu. \n Prob.", left_annotation = left_Annotation, 
                            bottom_annotation = bottom_Annotation, cluster_rows = F, 
                            cluster_columns = F, column_names_rot = 45, row_names_side = "left", 
                            col = col_map, split = as.character(group[sender.names]), 
                            column_split = as.character(group[receiver.names]), 
                            column_title = "Receiver", row_title = "Sender", 
                            column_title_side = "bottom", heatmap_legend_param = list(color_bar = "continuous"))
  }
  }  
 else{
   net_aggregated <- list()
  
   for (i in 1:length(dataset)) {
     print(dataset[i])
     if (is.null(sender.names)) {
       sender.names = rownames(object@net[[dataset[i]]][[1]])
     }
     if (is.null(receiver.names)) {
       receiver.names = colnames(object@net[[dataset[i]]][[1]])
     }
     if (interaction_use == "all") {
       net_aggregated[[i]] <- net_aggregation_v2(object@net[[dataset[i]]], method = method, 
                                            cut_off = cut_off)
     }
     else {
       net_aggregated[[i]] <- net_aggregation_v2(object@net[[dataset[i]]][interaction_use], 
                                            method = method, cut_off = cut_off)
     }
   }
   net_agg <- net_aggregated[[2]] - net_aggregated[[1]]
   # print("control")
   # print(net_aggregated[[1]][35:40,35:40])
   # print("MDD")
   # print(net_aggregated[[2]][35:40,35:40])
   # print("alternation")
   # print(net_agg[35:40,35:40])
   
   print("control")
   print(net_aggregated[[1]][36:41,36:41])
   print("MDD")
   print(net_aggregated[[2]][36:41,36:41])
   print("alternation")
   print(net_agg[36:41,36:41])
   # write.csv(net_aggregated[[1]], "mcontrol.csv")
   # write.csv(net_aggregated[[2]], "mmdd.csv")
   # write.csv(net_agg, "malter.csv")
   
   tmp_r <- rowSums(abs(net_agg))
   tmp_c <- colSums(abs(net_agg))
   tmp_df <- data.frame(
     cluster = names(tmp_r),
     sending = tmp_r,
     receiving = tmp_c,
     stringsAsFactors = FALSE
   )
   #write.csv(tmp_df[order(-tmp_df$receiving),],"rece.csv")
   #write.csv(tmp_df[order(-tmp_df$sending),],"send.csv")
   print(min(net_agg[sender.names, receiver.names]))
   print(max(net_agg[sender.names, receiver.names]))
   #print(rowSums(abs(net_agg)))
     #col_map = colorRamp2(c(min(net_agg[sender.names, 
      #                                            receiver.names]),0, max(net_agg[sender.names, 
       #                                                                                  receiver.names])), c("blue", "white", "red"))
     col_map = colorRamp2(c(-0.8,0, 0.3), c("#2461A1", "#F5F9FC", "#A91E2C")) # male: min -0.5, max 0.18; female: min -0.72, max 0.27
     if (is.null(group)) {
       ComplexHeatmap::Heatmap(net_agg[sender.names,
                                              receiver.names], name = method, cluster_rows = F,
                               cluster_columns = F, column_names_rot = 45, row_names_side = "left",
                               col = col_map, column_title = "Receiver", row_title = "Sender",
                               column_title_side = "bottom", heatmap_legend_param = list(color_bar = "continuous"))
     }
     else {
       left_Annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = scPalette(length(unique(group[sender.names])))),
                                                        labels = sort(unique(group[sender.names]))))
       bottom_Annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scPalette(length(unique(group[receiver.names])))),
                                                              labels = sort(unique(group[receiver.names]))))
       ComplexHeatmap::Heatmap(net_agg[sender.names,
                                              receiver.names], name = "Commu. \n Prob.", left_annotation = left_Annotation,
                               bottom_annotation = bottom_Annotation, cluster_rows = F,
                               cluster_columns = F, column_names_rot = 45, row_names_side = "left",
                               col = col_map, split = as.character(group[sender.names]),
                               column_split = as.character(group[receiver.names]),
                               column_title = "Receiver", row_title = "Sender",
                               column_title_side = "bottom", heatmap_legend_param = list(color_bar = "continuous"))
     }
 }
}