#' Plot functions for glsmr objects
#'
#' This function makes a few useful plots for glsmr objects
#'
#' @keywords plot glsmr
#' @param glsmr_obj a glsmr object from the function glsmr()
#' @param add_strata_2_curves binary TRUE or FALSE if strata boundaries should be illustrated to curve plots
#' @param add_strata_2_points binary TRUE or FALSE if strata boundaries should be illustrated to dot plots
#' @param brewer_col a brewer.pal color palette string
#' @return ggplot figure
#' @importFrom ggplot2 ggplot aes_string geom_point geom_rect geom_smooth aes scale_color_manual scale_shape_manual geom_errorbar geom_vline geom_hline labs theme_bw
#' @importFrom magrittr %>%
#' @importFrom ggpubr ggarrange
#' @importFrom dplyr mutate
#' @importFrom RColorBrewer brewer.pal
#' @export
#' @examples
#' plot_glsmr()
plot_glsmr = function( glsmr_obj, add_strata_2_curves = FALSE, add_strata_2_points = FALSE, brewer_col = "Set1" ){

  #############################
  ## I. Prepare SumStats for plots
  #############################
  ## plot data
  pdata = glsmr_obj$model_data

  ## table of summary statistics
  ss = glsmr_obj$summary_stats

  ## define variables
  outcome = ss$outcome
  exposure = ss$exposure
  instrument = "d.hat"
  strata_means = glsmr_obj$strata_IV_linear_mods$mean

  ## "full data" mean min and max for exposure
  fd_exp_mean = mean( glsmr_obj$model_data[, exposure], na.rm = TRUE)
  fd_exp_min = min( glsmr_obj$model_data[, exposure], na.rm = TRUE)
  fd_exp_max = max( glsmr_obj$model_data[, exposure], na.rm = TRUE)

  ## "full data" mean min and max for d.hat
  fd_dhat_mean = mean( glsmr_obj$model_data[, instrument], na.rm = TRUE)
  fd_dhat_min = min( glsmr_obj$model_data[, instrument], na.rm = TRUE)
  fd_dhat_max = max( glsmr_obj$model_data[, instrument], na.rm = TRUE)

  ## "full data" mean for outcome
  fd_outcome_mean = mean( glsmr_obj$model_data[, outcome], na.rm = TRUE)

  ## FORMAT P values for the linear vs non-linear F test
  obs_nonlinear_P = as.numeric( ss$obs_nonlinearity_test_P )
  obs_nonlinear_P = formatC(obs_nonlinear_P, format = "e", digits = 2)

  MR_nonlinear_P = as.numeric( ss$iv_nonlinearity_test_P )
  MR_nonlinear_P = formatC(MR_nonlinear_P, format = "e", digits = 2)

  ## Extract full observational model linear estimates
  l_est = summary(glsmr_obj$full_linear_model)$coef[exposure,]
  l_est = formatC(l_est, format = "e", digits = 2)

  ## Extract full ivreg (TSLS | MR) estimates
  iv_est = summary(glsmr_obj$full_ivreg_model )$coefficients[exposure,]
  iv_est = formatC(iv_est, format = "e", digits = 2)

  #############################
  ## II. Prepare data for strata
  ##      point estimate plots
  #############################
  ## linear estimate data
  ldata = glsmr_obj$strata_linear_mods[, c("beta","se","P","mean", "min","max")]
        ## defining the mean for min and max here as a geom_rect() plotting work around
  ldata = rbind(ldata, c(as.numeric(l_est[c(1,2,4)]), fd_exp_mean, fd_exp_mean, fd_exp_mean))
  ldata$data = "strata"
  ldata$data[nrow(ldata)] = "complete data set"
  ldata$data = as.factor(ldata$data)
  ldata$data = factor(ldata$data, levels = unique( ldata$data ) )

  ldata = ldata %>% mutate(sig = ifelse(P<0.05, "1", "0") )
  ldata$sig = factor(ldata$sig, levels = c("0","1"))

  ldata$strataID = rownames(ldata); ldata$strataID[nrow(ldata)] = "complete data set"
  ldata$strataID = factor(ldata$strataID, levels = ldata$strataID)

  ## ivreg estimate data
  ivdata = glsmr_obj$strata_ivreg_mods[, c("beta","se","P","mean", "min","max")]
        ## defining the mean for min and max here as a geom_rect() plotting work around
  ivdata = rbind(ivdata, c(as.numeric(iv_est[c(1,2,4)]), fd_exp_mean, fd_exp_mean, fd_exp_mean))
  ivdata$data = "strata"
  ivdata$data[nrow(ivdata)] = "complete data set"
  ivdata$data = as.factor(ivdata$data)
  ivdata$data = factor(ivdata$data, levels = unique( ivdata$data ) )

  ivdata = ivdata %>% mutate(sig = ifelse(P<0.05, "1", "0") )
  ivdata$sig = factor(ivdata$sig, levels = c("0","1"))

  ivdata$strataID = rownames(ivdata); ivdata$strataID[nrow(ivdata)] = "complete data set"
  ldata$strataID = factor(ldata$strataID, levels = ldata$strataID)

  ## linear IV estimate data: instrument-on-outcome / instrument-on-exposure
  l_iv_data = glsmr_obj$strata_IV_linear_mods[, c("beta_ratio","se_ratio","P","mean", "min","max")]
  colnames(l_iv_data) = c("beta","se","P","mean", "min","max")
      ## defining the mean for min and max here as a geom_rect() plotting work around
  l_iv_data = rbind(l_iv_data, c(as.numeric(iv_est[c(1,2,4)]), fd_exp_mean, fd_exp_mean, fd_exp_mean))
  l_iv_data$data = "strata"
  l_iv_data$data[nrow(l_iv_data)] = "complete data set"
  l_iv_data$data = as.factor(l_iv_data$data)
  l_iv_data$data = factor(l_iv_data$data, levels = unique( l_iv_data$data ) )

  l_iv_data = l_iv_data %>% mutate(sig = ifelse(P<0.05, "1", "0") )
  l_iv_data$sig = factor(l_iv_data$sig, levels = c("0","1"))

  l_iv_data$strataID = rownames(l_iv_data); l_iv_data$strataID[nrow(l_iv_data)] = "complete data set"
  l_iv_data$strataID = factor(l_iv_data$strataID, levels = l_iv_data$strataID)

  #############################
  ## III. PLOT linear non-linear
  ##     smooth plots
  #############################
  ## plot colors
  # pcol = c(RColorBrewer::brewer.pal(9, "Set1"), "grey20")

  ## plot linear and non-linear relationships
  o_lnl_plot = pdata %>% ggplot(aes_string(x = exposure, y = outcome)) +
    geom_smooth( method = "lm", formula = y~x, color = "black", se = TRUE ) +
    geom_smooth( method = "gam", formula = y~s(x), color = "blue", se = TRUE ) +
    geom_vline(xintercept = strata_means, color = "grey20", linetype = "dashed") +
    labs(title = "observational relationship",
         subtitle = paste0( "non-linear GAM a better fit? P = ", obs_nonlinear_P) ) +
    theme_bw()

  ## ADD strata color boundries
  if(add_strata_2_curves == TRUE){
    o_lnl_plot = o_lnl_plot +
      geom_rect(data = ldata[1:(nrow(ldata)-1),],
                aes(x = fd_exp_mean, y = fd_outcome_mean,
                    xmin = min, xmax = max,
                    ymin = -Inf, ymax = Inf,
                    fill = strataID ),
                alpha = 0.2 ) +
      scale_fill_brewer(palette = brewer_col)
  }

  iv_lnl_plot = pdata %>% ggplot(aes_string(x = instrument, y = outcome)) +
    geom_smooth( method = "lm", formula = y~x, color = "black", se = TRUE ) +
    geom_smooth( method = "gam", formula = y~s(x), color = "blue", se = TRUE ) +
    geom_vline(xintercept = strata_means, color = "grey20", linetype = "dashed") +
    labs(title = "MR relationship", x = paste0("genotype predicted ", exposure),
         subtitle = paste0( "non-linear GAM a better fit? P = ", MR_nonlinear_P) ) +
    theme_bw()

  ## ADD strata color boundries
  if(add_strata_2_curves == TRUE){
    iv_lnl_plot = iv_lnl_plot +
      geom_rect(data = ldata[1:(nrow(ldata)-1),],
                aes(x = fd_exp_mean, y = fd_outcome_mean,
                    xmin = min, xmax = max,
                    ymin = -Inf, ymax = Inf,
                    fill = strataID ),
                alpha = 0.2 ) +
      scale_fill_brewer(palette = brewer_col)
  }

  #############################
  ## IV. PLOT strata point
  ##     estimates
  #############################
  ## Strata point estimates
  if(add_strata_2_points == TRUE){
    l_est_plot = ldata %>%
      ggplot(aes(x = mean, y = beta)) +
      geom_rect(aes(xmin = min, xmax = max, ymin = -Inf, ymax = Inf, fill = strataID), alpha = 0.2) +
      scale_fill_brewer(palette = brewer_col) +
      geom_errorbar( aes(ymin=beta-se, ymax=beta+se), width = 0.2 ) +
      geom_point(aes(color = sig, shape = data ) , size = 3) +
      scale_color_manual(values = c("black","red"), drop = FALSE) +
      scale_shape_manual(values = c(19,22) ) +
      geom_hline(yintercept = 0, size = 1, color = "black", linetype = "dashed" ) +
      labs(x = "mean of exposure strata", title = "observational strata estimates") +
      theme_bw()

    iv_est_plot = l_iv_data %>%
      ggplot(aes(x = mean, y = beta)) +
      geom_rect(aes(xmin = min, xmax = max, ymin = -Inf, ymax = Inf, fill = strataID), alpha = 0.2) +
      scale_fill_brewer(palette = brewer_col) +
      geom_errorbar( aes(ymin=beta-se, ymax=beta+se), width = 0.2 ) +
      geom_point(aes(color = sig, shape = data  ) , size = 3) +
      scale_color_manual(values = c("black","red"), drop = FALSE) +
      scale_shape_manual(values = c(19,22) ) +
      geom_hline(yintercept = 0, size = 1, color = "black", linetype = "dashed" ) +
      labs(x = "mean of exposure strata", title = "tsls MR strata estimates") +
      theme_bw()
  } else {
    l_est_plot = ldata %>%
      ggplot(aes(x = mean, y = beta)) +
      geom_errorbar( aes(ymin=beta-se, ymax=beta+se), width = 0.2 ) +
      geom_point(aes(color = sig, shape = data ) , size = 3) +
      scale_color_manual(values = c("black","red"), drop = FALSE) +
      scale_shape_manual(values = c(19,22) ) +
      geom_hline(yintercept = 0, size = 1, color = "black", linetype = "dashed" ) +
      labs(x = "mean of exposure strata", title = "observational strata estimates") +
      theme_bw()

    iv_est_plot = l_iv_data %>%
      ggplot(aes(x = mean, y = beta)) +
      geom_errorbar( aes(ymin=beta-se, ymax=beta+se), width = 0.2 ) +
      geom_point(aes(color = sig, shape = data  ) , size = 3) +
      scale_color_manual(values = c("black","red"), drop = FALSE) +
      scale_shape_manual(values = c(19,22) ) +
      geom_hline(yintercept = 0, size = 1, color = "black", linetype = "dashed" ) +
      labs(x = "mean of exposure strata", title = "tsls MR strata estimates") +
      theme_bw()
  }

  #############################
  ## V. MERGE PLOTS
  #############################
  # plot1 = ggpubr::ggarrange(o_lnl_plot, iv_lnl_plot, nrow = 1 )
  plot1 = ggpubr::ggarrange(o_lnl_plot, iv_lnl_plot, nrow = 1, common.legend = TRUE,legend = "right" )
  plot2 = ggpubr::ggarrange(l_est_plot, iv_est_plot, nrow = 1, common.legend = TRUE,legend = "right" )

  plot = ggpubr::ggarrange(plot1, plot2, nrow = 2)

  #############################
  ## VI. Return to user
  #############################
  return(plot)

}
