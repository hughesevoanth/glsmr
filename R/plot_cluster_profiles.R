#' Plot the outcome exposure profiles derived from the function cluster_profiles()
#'
#' This function makes a plot of the outcome on exposure profiles
#'
#' @keywords GAM
#' @param cluster_data cluster_profiles() object
#' @param k which k (k-means) cluster to plot?
#' @param exposure_values what the x-axis values should be
#' @param exposure_name x-axis label
#' @param col_count how many columns to have in plot
#' @return a ggplot
#' @export
#' @examples
#' plot_cluster_profiles()
plot_cluster_profiles = function( cluster_data, k = 10,  exposure_values = NA , exposure_name = "Exposure", col_count = NULL){

  #############################################
  ## STEP 1:
  ##   Extract Data from cluster_profiles() object
  #############################################
  ## extract the clusters to plot
  myK = cluster_data$k[[k]]$cluster

  ## extract the predicted data matrix
  my_predicted_data = cluster_data$predictions

  #############################################
  ## STEP 2:
  ##   Define the average and SE of each profile
  ##   cluster
  #############################################
  ## Size of each cluster
  n = table(myK)
  ## estimate the mean and 95% CI of each cluster for each sample year month
  profiles = lapply(1:length(unique(myK)), function(j){
    w = which(myK == j)
    if(length(w)>1){
      ests = apply( my_predicted_data[w, ], 2, function(x){
        m = mean(x)
        s = quantile(x, probs = c(0.025,0.975))
        out = c(m, s)
        return(out)
      })
    } else {
      ests = rbind(my_predicted_data[w, ], my_predicted_data[w, ], my_predicted_data[w, ])
    }
    colnames(ests) = colnames(my_predicted_data)
    return(ests)
  })

  #############################################
  ## STEP 3:
  ##   Define the plot data frame
  #############################################
  ## Define what the (x-axis) exposure values should be
  ## if not defined
  if( is.na(exposure_values[1]) ){
    exposure_values = 1:ncol(my_predicted_data)
  }

  ## Define the plot data frame
  mytable = c()
  for(i in 1:length(profiles) ){
    x = t(profiles[[i]])
    colnames(x) = c("mean","lowerCI","upperCI")
    x = as.data.frame(x)
    x$cluster = paste0("cluster_",i, "_n=", n[i])
    x$exposure = exposure_values
    mytable = rbind(mytable, x)
  }

  ## Make sure the cluster column is a factor
  mytable$cluster = factor(mytable$cluster, levels = unique(mytable$cluster) )

  #############################################
  ## STEP 4:
  ##   Make the plot
  #############################################
  if(!is.null(col_count)){
    profile_plot = mytable %>% ggplot(aes( x = exposure, y = mean )) +
      geom_point(aes(color = as.factor(cluster) )) +
      geom_errorbar( aes(ymin= lowerCI , ymax= upperCI ),
                     width=.2,
                     position=position_dodge(0.05), color = "grey50" ) +
      facet_wrap(. ~ as.factor(cluster), ncol = col_count) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6)) +
      labs(color = "Kmeans Cluster", title = "Outcome on exposure GAM clustering profiles",
           y = "Outcome cluster means and 95% CI", x = exposure_name) +
      theme_bw() +
      theme(legend.position = "none")
  } else {
      profile_plot = mytable %>% ggplot(aes( x = exposure, y = mean )) +
        geom_point(aes(color = as.factor(cluster) )) +
        geom_errorbar( aes(ymin= lowerCI , ymax= upperCI ),
                       width=.2,
                       position=position_dodge(0.05), color = "grey50" ) +
        facet_wrap(. ~ as.factor(cluster)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6)) +
        labs(color = "Kmeans Cluster", title = "Outcome on exposure GAM clustering profiles",
             y = "Outcome cluster means and 95% CI", x = exposure_name) +
        theme_bw() +
        theme(legend.position = "none")
  }


  ## return to user
  return( profile_plot )


}
