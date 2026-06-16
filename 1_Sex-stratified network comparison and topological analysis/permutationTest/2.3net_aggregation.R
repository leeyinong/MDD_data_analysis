net_aggregation_v2 <- function (net_list, method = c("weight", "count", "weighted_count", 
                                                     "weighted_count2", "weight_threshold"), cut_off = 0.05) 
{
  method <- match.arg(method)
  if (method == "weighted_count") {
    net_aggregated <- 0 * net_list[[1]]
    for (jj in 1:length(net_list)) {
      net_aggregated <- net_aggregated + sum(net_list[[jj]]) * (net_list[[jj]] > 0)
    }
  }
  else if (method == "count") {
    net_aggregated <- 0 * net_list[[1]]
    for (jj in 1:length(net_list)) {
      net_aggregated <- net_aggregated + 1 * (net_list[[jj]] > 0)
    }
  }
  else if (method == "weighted_count2") {
    net_aggregated <- 0 * net_list[[1]]
    for (jj in 1:length(net_list)) {
      net_aggregated <- net_aggregated + sum(net_list[[jj]])/(1e-06 + sum(net_list[[jj]] > quantile(net_list[[jj]], probs = cut_off))) * (net_list[[jj]] > quantile(net_list[[jj]], probs = cut_off))
    }
  }
  else if (method == "weight_threshold") {
    net_aggregated <- 0 * net_list[[1]]
    for (jj in 1:length(net_list)) {
      net_aggregated <- net_aggregated + net_list[[jj]] * (net_list[[jj]] > quantile(net_list[[jj]], probs = cut_off))
    }
  }
  else {
    net_aggregated <- 0 * net_list[[1]]
    for (jj in 1:length(net_list)) {
      net_aggregated <- net_aggregated + net_list[[jj]]
    }
  }
  return(net_aggregated)
}