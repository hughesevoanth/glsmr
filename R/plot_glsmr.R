#' Plot functions for glsmr objects
#'
#' This function makes a few useful plots for glsmr objects
#'
#' @keywords plot glsmr
#' @param glsmr_obj a glsmr object from the function glsmr()
#' @return ggplot figure
#' @importFrom ggplot2 ggplot aes_string geom_point geom_smooth aes scale_color_manual scale_shape_manual geom_errorbar geom_vline geom_hline labs theme_bw
#' @importFrom magrittr %>%
#' @importFrom ggpubr ggarrange
#' @importFrom dplyr mutate
#' @export
#' @examples
#' plot_glsmr()
plot_glsmr = function( glsmr_obj ){

  #############################
  ## I. Prepare data for plots
  #############################
  ## plot data
  pdata = glsmr_obj$model_data

  ## table of summary statistics
  ss = as.data.frame(glsmr_obj$summary_stats)

  ## defive variables
  outcome = ss[ "outcome", ]
  exposure = ss[ "exposure", ]
  instrument = "d.hat"
  strata_means = glsmr_obj$strata_ivreg_mods$mean
  fulldata_mean = mean( glsmr_obj$model_data[, exposure], na.rm = TRUE)

  ## FORMAT P values for the linear vs non-linear F test
  obs_nonlinear_P = as.numeric( ss["linear_nonlinear_test_P",] )
  obs_nonlinear_P = formatC(obs_nonlinear_P, format = "e", digits = 2)

  MR_nonlinear_P = as.numeric( ss["iv_linear_nonlinear_test_P",] )
  MR_nonlinear_P = formatC(MR_nonlinear_P, format = "e", digits = 2)

  ## Extract full observational model linear estimates
  l_est = summary(glsmr_obj$full_linear_model)$coef[exposure,]
  l_est = formatC(l_est, format = "e", digits = 2)

  ## Extract full ivreg (TSLS | MR) estimates
  iv_est = summary(glsmr_obj$full_ivreg_model )$coefficients[exposure,]
  iv_est = formatC(iv_est, format = "e", digits = 2)

  #############################
  ## II. PLOT linear non-linear
  ##     smooth plots
  #############################
  ## plot linear and non-linear relationships
  o_lnl_plot = pdata %>% ggplot(aes_string(x = exposure, y = outcome)) +
    geom_smooth( method = "lm", formula = y~x, color = "black", se = TRUE ) +
    geom_smooth( method = "gam", formula = y~s(x), color = "blue", se = TRUE ) +
    geom_vline(xintercept = strata_means, color = "grey20", linetype = "dashed") +
    labs(title = "observational relationship",
         subtitle = paste0( "non-linear GAM a better fit? P = ", obs_nonlinear_P) ) +
    theme_bw()

  iv_lnl_plot = pdata %>% ggplot(aes_string(x = instrument, y = outcome)) +
    geom_smooth( method = "lm", formula = y~x, color = "black", se = TRUE ) +
    geom_smooth( method = "gam", formula = y~s(x), color = "blue", se = TRUE ) +
    geom_vline(xintercept = strata_means, color = "grey20", linetype = "dashed") +
    labs(title = "MR relationship", x = paste0("genotype predicted ", exposure),
         subtitle = paste0( "non-linear GAM a better fit? P = ", MR_nonlinear_P) ) +
    theme_bw()

  #############################
  ## III. Prepare data for strata
  ##      point estimate plots
  #############################
  ## linear estimate data
  ldata = glsmr_obj$strata_linear_mods[, c("beta","se","P","mean")]
  ldata = rbind(ldata, c(as.numeric(l_est[c(1,2,4)]), fulldata_mean))
  ldata$data = "strata"
  ldata$data[nrow(ldata)] = "all"
  ldata$data = as.factor(ldata$data)

  ldata = ldata %>% mutate(sig = ifelse(P<0.05, "1", "0") )
  ldata$sig = factor(ldata$sig, levels = c("0","1"))

  ## iv estimate data
  ivdata = glsmr_obj$strata_ivreg_mods[, c("beta","se","P","mean")]
  ivdata = rbind(ivdata, c(as.numeric(iv_est[c(1,2,4)]), fulldata_mean))
  ivdata$data = "strata"
  ivdata$data[nrow(ivdata)] = "all"
  ivdata$data = as.factor(ivdata$data)

  ivdata = ivdata %>% mutate(sig = ifelse(P<0.05, "1", "0") )
  ivdata$sig = factor(ivdata$sig, levels = c("0","1"))

  #############################
  ## IV. PLOT strata point
  ##     estimates
  #############################
  ## Strata point estimates
  l_est_plot = ldata %>%
    ggplot(aes(x = mean, y = beta)) +
    geom_errorbar( aes(ymin=beta-se, ymax=beta+se), width = 0.2 ) +
    geom_point(aes(color = sig, shape = data ) , size = 3) +
    scale_color_manual(values = c("black","red"), drop = FALSE) +
    scale_shape_manual(values = c(22,19) ) +
    geom_hline(yintercept = 0, size = 1, color = "black", linetype = "dashed" ) +
    labs(x = "mean of exposure strata", title = "observational strata estimates") +
    theme_bw()

  iv_est_plot = ivdata %>%
    ggplot(aes(x = mean, y = beta)) +
    geom_errorbar( aes(ymin=beta-se, ymax=beta+se), width = 0.2 ) +
    geom_point(aes(color = sig, shape = data  ) , size = 3) +
    scale_color_manual(values = c("black","red"), drop = FALSE) +
    scale_shape_manual(values = c(22,19) ) +
    geom_hline(yintercept = 0, size = 1, color = "black", linetype = "dashed" ) +
    labs(x = "mean of exposure strata", title = "tsls MR strata estimates") +
    theme_bw()

  #############################
  ## V. MERGE PLOTS
  #############################
  plot1 = ggpubr::ggarrange(o_lnl_plot, iv_lnl_plot, nrow = 1 )
  plot2 = ggpubr::ggarrange(l_est_plot, iv_est_plot, nrow = 1, common.legend = TRUE,legend = "right" )

  plot = ggpubr::ggarrange(plot1, plot2, nrow = 2)

  #############################
  ## VI. Return to user
  #############################
  return(plot)

}
