#' K-means Clustering of OUtcome on Exposure profiles
#'
#' This function takes a matrix of outcome on exposure predicted values and performs a kmeans clustering to aid in identifying similar profiles
#'
#' @keywords GAM
#' @param wdata a data frame of predicted exposure values, such as those derived from the function gam_predict()
#' @return a list of kmeans values, and plots
#' @importFrom stats kmeans
#' @importFrom ggpubr ggarrange
#' @importFrom cluster silhouette
#' @export
#' @examples
#' cluster_profiles()
cluster_profiles = function( wdata,
                             center = TRUE,
                             scale = TRUE,
                             k_count = 20
                  ){
  #############################################
  ## STEP 1:
  ## mean center and scale the predictions
  #############################################
  predictions = t( apply(wdata, 1, function(x){ scale(x, center = center, scale = scale) }) )
  colnames(predictions) = colnames(wdata)

  #############################################
  ## Step 2:
  ## run the kmeans analysis
  #############################################
  set.seed(1234)
  k = lapply(2:k_count, function(k){
    out = kmeans( predictions , centers = k, iter.max = 500, nstart = 50)
    return(out)
  } )

  #############################################
  ## Step 3:
  ##  Run Silhouette, extract widths and betweeness and withness
  #############################################
  ## Distance matrix
  dmat = dist( predictions )
  ## extract widths
  silhouette_score = unlist( lapply(k, function(x){
    ss = cluster::silhouette(x$cluster, dmat)
    score = mean(ss[, 3])
  }) )
  ## Extract BETWEEN and Withiness Sums of Squares
  BSS = unlist( lapply(k, function(x){ x$betweenss / x$totss}) )
  WSS = unlist( lapply(k, function(x){ x$tot.withinss / x$totss}) )

  #############################################
  ## Step 3:
  ##  Generate Plots
  #############################################
  ### Define Plot Data
  plotdata = data.frame(k = 1:k_count, swidth = c(0, silhouette_score), BSS = c(0,BSS), WSS = c(0, WSS) )
  ####
  plot1 = plotdata[-1,] %>% ggplot(aes(x = k, y = swidth)) +
    geom_point( color = "blue", size = 3) +
    geom_line() +
    geom_text(aes(label = k), vjust = -1, size = 3) +
    geom_vline(xintercept = which(plotdata[,2] == max(plotdata[,2]) ) , color = "red" ) +
    labs(x = 'Number of Clusters', y = 'Average Silhouette Scores') +
    theme_bw()

  plot2 = plotdata[-1, ] %>% ggplot(aes(x = k, y = WSS) ) +
    geom_point( color = "blue", size = 3) +
    geom_line() +
    geom_text(aes(label = k), vjust = -1, size = 3) +
    geom_vline(xintercept = which(plotdata[,2] == max(plotdata[,2]) ) , color = "red" ) +
    labs(x = 'Number of Clusters', y = 'Within Cluster Sum of Squares / Total SS') +
    theme_bw()

  plots = ggpubr::ggarrange(plot1, plot2, nrow = 1)

  #############################################
  ## Step 4:
  ##  Return Results
  #############################################
  out = list(k = k,
             plots = plots,
             predictions = predictions )
  return(out)


}
