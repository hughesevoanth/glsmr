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
#' @importFrom ggplot2 ggplot aes_string geom_point geom_rect geom_smooth aes scale_color_manual scale_shape_manual scale_fill_brewer geom_errorbar geom_vline geom_hline labs theme_bw
#' @importFrom magrittr %>%
#' @importFrom ggpubr ggarrange
#' @importFrom dplyr mutate
#' @export
#' @examples
#' plot_glsmr()
plot_glsmr = function( glsmr_obj, add_strata_2_curves = FALSE, add_strata_2_points = FALSE, brewer_col = "Set1" ){

  #############################
  ## I. Prepare SumStats for plots
  #############################
  ## plot data
  pdata = na.omit(glsmr_obj$summary_tables$model_data)

  ## table of summary statistics
  ss = glsmr_obj$summary_stats

  ## define variables
  outcome = as.character( ss$outcome )
  exposure = as.character(ss$exposure )
  instrument = as.character( ss$instrument )
  iv = as.character( "iv_predicted_exposure" )

  ## summary statistics for strata and full data models
  obs_strata_ss = glsmr_obj$summary_tables$observational_sumstats
  tsls_strata_ratio_ss = glsmr_obj$summary_tables$strata_tsls_ratio_sumstats
  tsls_strata_ivreg_ss = glsmr_obj$summary_tables$strata_tsls_ivreg_sumstats

  #############################
  ## II. Prepare data for strata
  ##      point estimate plots
  #############################
  ## Observational table
  obs_strata_ss$sig = ifelse(obs_strata_ss$P <= 0.05, "1", "0"  )
  obs_strata_ss$sig = factor(obs_strata_ss$sig, levels = c("0","1"))
  obs_strata_ss$strataID = rownames(obs_strata_ss); obs_strata_ss$strataID[nrow(obs_strata_ss)] = "complete data set"
  obs_strata_ss$strataID = factor(obs_strata_ss$strataID, levels = obs_strata_ss$strataID)
  obs_strata_ss$shape = "strata"; obs_strata_ss$shape[nrow(obs_strata_ss)] = "complete data set"
  obs_strata_ss$shape = factor(obs_strata_ss$shape, levels = c("strata","complete data set") )
  obs_strata_ss["fulldata","min"] = obs_strata_ss["fulldata","mean"]
  obs_strata_ss["fulldata","max"] = obs_strata_ss["fulldata","mean"]
  ## RATIO TABLE
  tsls_strata_ratio_ss$sig = ifelse(tsls_strata_ratio_ss$P <= 0.05, "1", "0"  )
  tsls_strata_ratio_ss$sig = factor(tsls_strata_ratio_ss$sig, levels = c("0","1"))
  tsls_strata_ratio_ss$strataID = rownames(tsls_strata_ratio_ss); tsls_strata_ratio_ss$strataID[nrow(tsls_strata_ratio_ss)] = "complete data set"
  tsls_strata_ratio_ss$strataID = factor(tsls_strata_ratio_ss$strataID, levels = tsls_strata_ratio_ss$strataID)
  tsls_strata_ratio_ss$shape = "strata"; tsls_strata_ratio_ss$shape[nrow(tsls_strata_ratio_ss)] = "complete data set"
  tsls_strata_ratio_ss$shape = factor(tsls_strata_ratio_ss$shape, levels = c("strata","complete data set") )
  tsls_strata_ratio_ss["fulldata","min"] = tsls_strata_ratio_ss["fulldata","mean"]
  tsls_strata_ratio_ss["fulldata","max"] = tsls_strata_ratio_ss["fulldata","mean"]
  ## IVREG TABLE
  tsls_strata_ivreg_ss$sig = ifelse(tsls_strata_ivreg_ss$P <= 0.05, "1", "0"  )
  tsls_strata_ivreg_ss$sig = factor(tsls_strata_ivreg_ss$sig, levels = c("0","1"))
  tsls_strata_ivreg_ss$strataID = rownames(tsls_strata_ivreg_ss); tsls_strata_ivreg_ss$strataID[nrow(tsls_strata_ivreg_ss)] = "complete data set"
  tsls_strata_ivreg_ss$strataID = factor(tsls_strata_ivreg_ss$strataID, levels = tsls_strata_ivreg_ss$strataID)
  tsls_strata_ivreg_ss$shape = "strata"; tsls_strata_ivreg_ss$shape[nrow(tsls_strata_ivreg_ss)] = "complete data set"
  tsls_strata_ivreg_ss$shape = factor(tsls_strata_ivreg_ss$shape, levels = c("strata","complete data set") )
  tsls_strata_ivreg_ss["fulldata","min"] = tsls_strata_ivreg_ss["fulldata","mean"]
  tsls_strata_ivreg_ss["fulldata","max"] = tsls_strata_ivreg_ss["fulldata","mean"]

  ## "full data" mean min and max for exposure
  fd_exp_mean = mean( pdata[, exposure], na.rm = TRUE)

  ## "full data" mean min and max for d.hat
  fd_dhat_mean = mean( pdata[, iv], na.rm = TRUE)

  ## "full data" mean for outcome
  fd_outcome_mean = mean( pdata[, outcome], na.rm = TRUE)

  ## FORMAT P values for the linear vs non-linear F test
  obs_nonlinear_P = formatC( glsmr_obj$summary_tables$nonlinearity_test["obs", "P"] , format = "e", digits = 2)
  MR_nonlinear_P = formatC( glsmr_obj$summary_tables$nonlinearity_test["tsls", "P"] , format = "e", digits = 2)

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
    geom_point( x = fd_exp_mean , y = fd_outcome_mean , size = 3, color = "darkorange1", shape = 15 ) +
    geom_vline(xintercept = obs_strata_ss$mean[ -length(obs_strata_ss$mean) ], color = "grey20", linetype = "dashed") +
    labs(title = "observational relationship",
         subtitle = paste0( "non-linear GAM a better fit? P = ", obs_nonlinear_P) ) +
    theme_bw()

  ## ADD strata color boundries
  if(add_strata_2_curves == TRUE){
    o_lnl_plot = o_lnl_plot +
      geom_rect(data = obs_strata_ss[1:(nrow(obs_strata_ss)-1),],
                aes(x = NULL, y = NULL,
                    xmin = min, xmax = max,
                    ymin = -Inf, ymax = Inf,
                    fill = strataID ),
                alpha = 0.2 ) +
      scale_fill_brewer(palette = brewer_col)
  }

  iv_lnl_plot = pdata %>% ggplot(aes_string(x = iv, y = outcome)) +
    geom_smooth( method = "lm", formula = y~x, color = "black", se = TRUE ) +
    geom_smooth( method = "gam", formula = y~s(x), color = "blue", se = TRUE ) +
    geom_point( aes(x = fd_dhat_mean , y = fd_outcome_mean) , size = 3, color = "darkorange1", shape = 15 ) +
    # geom_vline(xintercept = tsls_strata_ratio_ss$mean[ -length(tsls_strata_ratio_ss$mean) ], color = "grey20", linetype = "dashed") +
    labs(title = "MR relationship", x = paste0("genotype predicted ", exposure),
         subtitle = paste0( "non-linear GAM a better fit? P = ", MR_nonlinear_P) ) +
    theme_bw()

  ## ADD strata color boundries
  # if(add_strata_2_curves == TRUE){
  #   iv_lnl_plot = iv_lnl_plot +
  #     geom_rect(data = tsls_strata_ratio_ss[ 1:(nrow(tsls_strata_ratio_ss)-1) , ],
  #               aes(x = NULL, y = NULL,
  #                   xmin = min, xmax = max,
  #                   ymin = -Inf, ymax = Inf,
  #                   fill = strataID ),
  #               alpha = 0.2 ) +
  #     scale_fill_brewer(palette = brewer_col)
  # }

  #############################
  ## IV. PLOT strata point
  ##     estimates
  #############################
  ## Strata point estimates
  if(add_strata_2_points == TRUE){
    l_est_plot = obs_strata_ss %>%
      ggplot(aes(x = mean, y = beta)) +
      geom_rect(data = obs_strata_ss[1:(nrow(tsls_strata_ratio_ss)-1),],
                aes(xmin = min, xmax = max, ymin = -Inf, ymax = Inf, fill = strataID),
                alpha = 0.2) +
      scale_fill_brewer(palette = brewer_col) +
      geom_errorbar( aes(ymin=beta-(1.96*se), ymax=beta+(1.96*se)), width = 0.2 ) +
      geom_point(aes(color = sig, shape = shape ) , size = 3) +
      scale_color_manual(values = c("black","red"), drop = FALSE) +
      scale_shape_manual(values = c(19,22) ) +
      geom_hline(yintercept = 0, size = 1, color = "black", linetype = "dashed" ) +
      labs(x = "mean of exposure strata", title = "observational strata estimates") +
      theme_bw()

    iv_est_plot = tsls_strata_ratio_ss %>%
      ggplot(aes(x = mean, y = beta_iv)) +
      geom_rect(data = tsls_strata_ratio_ss[1:(nrow(tsls_strata_ratio_ss)-1),],
                aes(xmin = min, xmax = max, ymin = -Inf, ymax = Inf, fill = strataID),
                alpha = 0.2) +
      scale_fill_brewer(palette = brewer_col) +
      geom_errorbar( aes(ymin=beta_iv-(1.96*se_iv), ymax=beta_iv+(1.96*se_iv)), width = 0.2 ) +
      geom_point(aes(color = sig, shape = shape  ) , size = 3) +
      scale_color_manual(values = c("black","red"), drop = FALSE) +
      scale_shape_manual(values = c(19,22) ) +
      geom_hline(yintercept = 0, size = 1, color = "black", linetype = "dashed" ) +
      labs(x = "mean of exposure strata", y = "beta" ,title = "tsls MR strata estimates") +
      theme_bw()
  } else {
    l_est_plot = obs_strata_ss %>%
      ggplot(aes(x = mean, y = beta)) +
      geom_errorbar( aes(ymin=beta-(1.96*se), ymax=beta+(1.96*se)), width = 0.2 ) +
      geom_point(aes(color = sig, shape = shape ) , size = 3) +
      scale_color_manual(values = c("black","red"), drop = FALSE) +
      scale_shape_manual(values = c(19,22) ) +
      geom_hline(yintercept = 0, size = 1, color = "black", linetype = "dashed" ) +
      labs(x = "mean of exposure strata", title = "observational strata estimates") +
      theme_bw()

    iv_est_plot = tsls_strata_ratio_ss %>%
      ggplot(aes(x = mean, y = beta_iv)) +
      geom_errorbar( aes(ymin=beta_iv-(1.96*se_iv), ymax=beta_iv+(1.96*se_iv)), width = 0.2 ) +
      geom_point(aes(color = sig, shape = shape  ) , size = 3) +
      scale_color_manual(values = c("black","red"), drop = FALSE) +
      scale_shape_manual(values = c(19,22) ) +
      geom_hline(yintercept = 0, size = 1, color = "black", linetype = "dashed" ) +
      labs(x = "mean of exposure strata", y = "beta", title = "tsls MR strata estimates") +
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
