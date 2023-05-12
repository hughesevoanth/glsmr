#' Plot functions for glsmr objects
#'
#' This function makes a few useful plots for glsmr objects
#'
#' @keywords plot glsmr
#' @param glsmr_obj a glsmr object from the function glsmr()
#' @param add_strata_2_curves binary TRUE or FALSE if strata boundaries should be illustrated to curve plots
#' @param add_strata_2_points binary TRUE or FALSE if strata boundaries should be illustrated to dot plots
#' @param brewer_col a brewer.pal color palette string
#' @param old_plot_scheme set to TRUE to get an older version of the plot
#' @param old_GAM_smooths binary TRUE or FALSE. The Default is FALSE. If TRUE the linear and GAM smooths will be derived from the complete data set and the geom_smooth() function of ggplot2, not from predicted values derived from the fitted models.
#' @param plot_obs_res_betas set to TRUE to plot GAM residualized effect estimates for observational data
#' @param pval_thresh The pvalue threshold to use to indicate which point estimates indicate an association between exposure and outcome. The default value is 0.05.
#' @return ggplot figure
#' @importFrom ggplot2 ggplot aes_string geom_point geom_rect geom_smooth aes scale_color_manual scale_shape_manual scale_fill_brewer geom_errorbar geom_vline geom_hline labs theme_bw
#' @importFrom magrittr %>%
#' @importFrom ggpubr ggarrange
#' @importFrom dplyr mutate
#' @importFrom ggnewscale new_scale
#' @export
#' @examples
#' plot_glsmr()
plot_glsmr = function( glsmr_obj,
                       add_strata_2_curves = FALSE,
                       add_strata_2_points = FALSE,
                       brewer_col = "Set1",
                       old_plot_scheme = FALSE,
                       old_GAM_smooths = FALSE,
                       plot_obs_res_betas = FALSE,
                       pval_thresh = 0.05){

  #############################
  ## I. Prepare SumStats for plots
  #############################
  ## plot data
  pdata = na.omit(glsmr_obj$summary_tables$model_data)

  ## table of summary statistics
  ss = glsmr_obj$summary_stats

  ## non-linearity tests
  nl_tests = glsmr_obj$summary_tables$nonlinearity_test

  ## meta tests
  meta_test = glsmr_obj$summary_tables$strata_meta_analysis

  ## define variables
  outcome = as.character( ss$outcome )
  exposure = as.character(ss$exposure )
  instrument = as.character( ss$instrument )
  iv = as.character( "iv_predicted_exposure" )

  ## summary statistics for strata and full data models
  obs_strata_ss = glsmr_obj$summary_tables$observational_sumstats
  obs_res_strata_ss = glsmr_obj$summary_tables$observational_GAM_res_sumstats
  tsls_strata_ratio_ss = glsmr_obj$summary_tables$strata_tsls_ratio_sumstats
  tsls_strata_ivreg_ss = glsmr_obj$summary_tables$strata_tsls_ivreg_sumstats

  #############################
  ## II. Prepare data for strata
  ##      point estimate plots
  #############################
  ## define pvalue string
  if( nchar(pval_thresh) > 4 ){
    pv = formatC(pval_thresh, format = "e", digits = 2)
  } else {
    pv = as.character( pval_thresh )
    }
  ## Observational table
  obs_strata_ss$sig = ifelse(obs_strata_ss$P <= pval_thresh, paste0("P <= ", pv ) ,  paste0("P > ", pv )  )
  obs_strata_ss$sig = factor(obs_strata_ss$sig, levels = c( paste0("P > ", pv ) ,  paste0("P <= ", pv ) ) )
  obs_strata_ss$strataID = rownames(obs_strata_ss); obs_strata_ss$strataID[nrow(obs_strata_ss)] = "complete data set"
  obs_strata_ss$strataID = factor(obs_strata_ss$strataID, levels = obs_strata_ss$strataID)
  obs_strata_ss$shape = "strata"; obs_strata_ss$shape[nrow(obs_strata_ss)] = "complete data set"
  obs_strata_ss$shape = factor(obs_strata_ss$shape, levels = c("strata","complete data set") )
  obs_strata_ss["fulldata","exposure_min"] = obs_strata_ss["fulldata","exposure_mean"]
  obs_strata_ss["fulldata","exposure_max"] = obs_strata_ss["fulldata","exposure_mean"]

  ## Residual Observational table
  obs_res_strata_ss$sig = ifelse(obs_res_strata_ss$P <= pval_thresh, paste0("P <= ", pv ) ,  paste0("P > ", pv )  )
  obs_res_strata_ss$sig = factor(obs_res_strata_ss$sig, levels = c( paste0("P > ", pv ) ,  paste0("P <= ", pv ) ) )
  obs_res_strata_ss$strataID = rownames(obs_res_strata_ss); obs_res_strata_ss$strataID[nrow(obs_res_strata_ss)] = "complete data set"
  obs_res_strata_ss$strataID = factor(obs_res_strata_ss$strataID, levels = obs_res_strata_ss$strataID)
  obs_res_strata_ss$shape = "strata"; obs_res_strata_ss$shape[nrow(obs_res_strata_ss)] = "complete data set"
  obs_res_strata_ss$shape = factor(obs_res_strata_ss$shape, levels = c("strata","complete data set") )
  obs_res_strata_ss["fulldata","exposure_min"] = obs_res_strata_ss["fulldata","exposure_mean"]
  obs_res_strata_ss["fulldata","exposure_max"] = obs_res_strata_ss["fulldata","exposure_mean"]

  ## RATIO TABLE
  tsls_strata_ratio_ss$sig = ifelse(tsls_strata_ratio_ss$P <= pval_thresh, paste0("P <= ", pv ) ,  paste0("P > ", pv ) )
  tsls_strata_ratio_ss$sig = factor(tsls_strata_ratio_ss$sig, levels = c( paste0("P > ", pv ) ,  paste0("P <= ", pv ) ) )
  tsls_strata_ratio_ss$strataID = rownames(tsls_strata_ratio_ss); tsls_strata_ratio_ss$strataID[nrow(tsls_strata_ratio_ss)] = "complete data set"
  tsls_strata_ratio_ss$strataID = factor(tsls_strata_ratio_ss$strataID, levels = tsls_strata_ratio_ss$strataID)
  tsls_strata_ratio_ss$shape = "strata"; tsls_strata_ratio_ss$shape[nrow(tsls_strata_ratio_ss)] = "complete data set"
  tsls_strata_ratio_ss$shape = factor(tsls_strata_ratio_ss$shape, levels = c("strata","complete data set") )
  tsls_strata_ratio_ss["fulldata","exposure_min"] = tsls_strata_ratio_ss["fulldata","exposure_mean"]
  tsls_strata_ratio_ss["fulldata","exposure_max"] = tsls_strata_ratio_ss["fulldata","exposure_mean"]

  ## IVREG TABLE
  tsls_strata_ivreg_ss$sig = ifelse(tsls_strata_ivreg_ss$P <= pval_thresh, paste0("P <= ", pv ) ,  paste0("P > ", pv ) )
  tsls_strata_ivreg_ss$sig = factor(tsls_strata_ivreg_ss$sig, levels = c( paste0("P > ", pv ) ,  paste0("P <= ", pv ) ) )
  tsls_strata_ivreg_ss$strataID = rownames(tsls_strata_ivreg_ss); tsls_strata_ivreg_ss$strataID[nrow(tsls_strata_ivreg_ss)] = "complete data set"
  tsls_strata_ivreg_ss$strataID = factor(tsls_strata_ivreg_ss$strataID, levels = tsls_strata_ivreg_ss$strataID)
  tsls_strata_ivreg_ss$shape = "strata"; tsls_strata_ivreg_ss$shape[nrow(tsls_strata_ivreg_ss)] = "complete data set"
  tsls_strata_ivreg_ss$shape = factor(tsls_strata_ivreg_ss$shape, levels = c("strata","complete data set") )
  tsls_strata_ivreg_ss["fulldata","exposure_min"] = tsls_strata_ivreg_ss["fulldata","exposure_mean"]
  tsls_strata_ivreg_ss["fulldata","exposure_max"] = tsls_strata_ivreg_ss["fulldata","exposure_mean"]

  ## "full data" mean min and max for exposure
  fd_exp_mean = mean( pdata[, exposure], na.rm = TRUE)

  ## "full data" mean min and max for d.hat
  fd_dhat_mean = mean( pdata[, iv], na.rm = TRUE)

  ## "full data" mean for outcome
  fd_outcome_mean = mean( pdata[, outcome], na.rm = TRUE)

  ## FORMAT P values for the linear vs non-linear F test
  obs_nonlinear_P = formatC( nl_tests["obs", "obs_nonlinear_Ftest_P"] , format = "e", digits = 2)
  MR_nonlinear_P = formatC( nl_tests["tsls", "obs_nonlinear_Ftest_P"] , format = "e", digits = 2)

  ## LRT
  temp_p = ifelse(nl_tests["obs", "obs_nonlinear_LRT_ChisqStat"] < 0.5, 1, nl_tests["obs", "obs_nonlinear_LRT_P"])
  obs_nonlinear_LRT_P = formatC( temp_p , format = "e", digits = 2)
  temp_p = ifelse(nl_tests["tsls", "obs_nonlinear_LRT_ChisqStat"] < 0.5, 1, nl_tests["tsls", "obs_nonlinear_LRT_P"])
  MR_nonlinear_LRT_P = formatC( temp_p , format = "e", digits = 2)

  ## FORMAT P values for the heterogenity tests
  obs_het_P = formatC( meta_test["obs", "meta_het_p"] , format = "e", digits = 2)
  res_obs_het_P = formatC( meta_test["res_obs", "meta_het_p"] , format = "e", digits = 2)
  MR_het_P = formatC( meta_test["mr", "meta_het_p"] , format = "e", digits = 2)

  ## FORMAT P values for the trend tests
  obs_trend_P = formatC( meta_test["obs", "meta_strata_mean_p"] , format = "e", digits = 2)
  res_obs_trend_P = formatC( meta_test["res_obs", "meta_strata_mean_p"] , format = "e", digits = 2)
  MR_trend_P = formatC( meta_test["mr", "meta_strata_mean_p"] , format = "e", digits = 2)


  #############################
  ## III. Extract Predictions
  ##     from GAM Models
  #############################

  ## (A) OBSERVATIONAL PREDICTIONS ****************

  ## Extract Predictions for the NULL (exposure is modeled as a parametric term) GAM model.
  predicted_tbl = tidymv::get_gam_predictions(glsmr_obj$models_no_stratification$observational_gam_null,
                                  series = {{exposure}},
                                  series_length = 45 )

  ## identify the other covariates in the predicted table
  ovars = colnames(predicted_tbl)
  w = which(!ovars %in% c( ".idx", "SE", "CI_upper", "CI_lower", exposure, outcome) )
  if(length(w)>0){
    linear_covariates = ovars[w]
  } else { linear_covariates = NA }

  ## limit to one value of each linear covariate
  if(!is.na(linear_covariates[1])){
    for(lc in linear_covariates){
      w = which(colnames(predicted_tbl) == lc)
      vals = unique(predicted_tbl[,w])
      if(length(vals) == 2){ ## if length is 2 as it might be for sex, just take the first one.
        vals = vals[1]
      } else { ## take the middle value
        vals = vals[ round( mean( 1:length(vals) ) ) ]
      }
      ## Reduce
      predicted_tbl = predicted_tbl %>% filter( .data[[lc]] == vals )
    }
  }

  ## Add model type to the table
  predicted_tbl$model = "parametric"
  predicted_tbl_P = predicted_tbl
  rm(predicted_tbl)

  ## Extract Predictions for the GAM (exposure is modeled as a smooth term) model.
  predicted_tbl = tidymv::get_gam_predictions(glsmr_obj$models_no_stratification$observational_gam,
                                  series = {{exposure}} ,
                                  series_length = 45 )

  ## identify the other covariates in the predicted table
  ovars = colnames(predicted_tbl)
  w = which(!ovars %in% c( ".idx", "SE", "CI_upper", "CI_lower", exposure, outcome) )
  if(length(w)>0){
    linear_covariates = ovars[w]
  } else { linear_covariates = NA }

  ## limit to one value of each linear covariate
  if(!is.na(linear_covariates[1])){
    for(lc in linear_covariates){
      w = which(colnames(predicted_tbl) == lc)
      vals = unique(predicted_tbl[,w])
      if(length(vals) == 2){ ## if length is 2 as it might be for sex, just take the first one.
        vals = vals[1]
      } else { ## take the middle value
        vals = vals[ round( mean( 1:length(vals) ) ) ]
      }
      ## Reduce
      predicted_tbl = predicted_tbl %>% filter( .data[[lc]] == vals )
    }
  }
  ## Add model type to the table
  predicted_tbl$model = "smooth"
  predicted_tbl_S = predicted_tbl
  rm(predicted_tbl)

  ## Pull the Parametric and Smooth Data together
  obs_predicted_tbl = rbind(predicted_tbl_P, predicted_tbl_S)


  ## (B) MR PREDICTIONS ****************

  ## Extract Predictions for the NULL (exposure is modeled as a parametric term) GAM model.
  predicted_tbl = tidymv::get_gam_predictions(glsmr_obj$models_no_stratification$tsls_gam_null,
                                              series = {{iv}},
                                              series_length = 45 )

  ## identify the other covariates in the predicted table
  ovars = colnames(predicted_tbl)
  w = which(!ovars %in% c( ".idx", "SE", "CI_upper", "CI_lower", iv, outcome) )
  if(length(w)>0){
    linear_covariates = ovars[w]
  } else { linear_covariates = NA }

  ## limit to one value of each linear covariate
  if(!is.na(linear_covariates[1])){
    for(lc in linear_covariates){
      w = which(colnames(predicted_tbl) == lc)
      vals = unique(predicted_tbl[,w])
      if(length(vals) == 2){ ## if length is 2 as it might be for sex, just take the first one.
        vals = vals[1]
      } else { ## take the middle value
        vals = vals[ round( mean( 1:length(vals) ) ) ]
      }
      ## Reduce
      predicted_tbl = predicted_tbl %>% filter( .data[[lc]] == vals )
    }
  }
  ## Add model type to the table
  predicted_tbl$model = "parametric"
  predicted_tbl_P = predicted_tbl
  rm(predicted_tbl)

  ## Extract Predictions for the GAM (exposure is modeled as a smooth term) model.
  predicted_tbl = tidymv::get_gam_predictions(glsmr_obj$models_no_stratification$tsls_gam,
                                              series = {{iv}} ,
                                              series_length = 45 )

  ## identify the other covariates in the predicted table
  ovars = colnames(predicted_tbl)
  w = which(!ovars %in% c( ".idx", "SE", "CI_upper", "CI_lower", iv, outcome) )
  if(length(w)>0){
    linear_covariates = ovars[w]
  } else { linear_covariates = NA }

  ## limit to one value of each linear covariate
  if(!is.na(linear_covariates[1])){
    for(lc in linear_covariates){
      w = which(colnames(predicted_tbl) == lc)
      vals = unique(predicted_tbl[,w])
      if(length(vals) == 2){ ## if length is 2 as it might be for sex, just take the first one.
        vals = vals[1]
      } else { ## take the middle value
        vals = vals[ round( mean( 1:length(vals) ) ) ]
      }
      ## Reduce
      predicted_tbl = predicted_tbl %>% filter( .data[[lc]] == vals )
    }
  }
  ## Add model type to the table
  predicted_tbl$model = "smooth"
  predicted_tbl_S = predicted_tbl
  rm(predicted_tbl)

  ## Pull the Parametric and Smooth Data together
  iv_predicted_tbl = rbind(predicted_tbl_P, predicted_tbl_S)


  #############################
  ## IV. PLOT linear non-linear
  ##     smooth plots
  #############################
  if(old_GAM_smooths == TRUE){
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
      o_lnl_plot = o_lnl_plot + new_scale("fill") +
        geom_rect(data = obs_strata_ss[1:(nrow(obs_strata_ss)-1),],
                  aes(x = NULL, y = NULL,
                      xmin = exposure_min, xmax = exposure_max,
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

  } else {
    ## NEW PLOT SCHEME WITH GAM PREDICTION CURVES
    ## OBSERVATIONAL PLOT
    o_lnl_plot = obs_predicted_tbl %>% ggplot( aes_string(x = exposure, y = outcome) ) +
      geom_line( aes( color = model ), size = 1.5) +
      geom_ribbon(aes(ymin = CI_lower,  ymax = CI_upper, fill = model ), alpha = 0.2 ) +
      #scale_color_brewer(palette = "Set1") +
      #scale_fill_brewer(palette = "Set1") +
      scale_color_manual(values = c("black", "blue")) +
      scale_fill_manual(values = c("black", "blue")) +
      geom_point( x = fd_exp_mean , y = fd_outcome_mean , size = 3, color = "darkorange1", shape = 15 ) +
      labs(title = "observational relationship",
           subtitle = paste0( "parametic-v-smooth: F-test P = ", obs_nonlinear_P, "; LRT P = ", obs_nonlinear_LRT_P) ) +
      theme_bw() +
      theme(legend.position = "none")

    ## ADD strata color boundries
    if(add_strata_2_curves == TRUE){
      o_lnl_plot = o_lnl_plot + new_scale("fill") +
        geom_rect(data = obs_strata_ss[1:(nrow(obs_strata_ss)-1),],
                  aes(x = NULL, y = NULL,
                      xmin = exposure_min, xmax = exposure_max,
                      ymin = -Inf, ymax = Inf,
                      fill = strataID ),
                  alpha = 0.2 ) +
        scale_fill_brewer(palette = brewer_col) +
        theme_bw()
    }

    ## MR | TSLS PLOT
    iv_lnl_plot = iv_predicted_tbl %>% ggplot( aes_string(x = iv, y = outcome) ) +
      geom_line( aes( color = model ) , size = 1.5) +
      geom_ribbon(aes(ymin = CI_lower,  ymax = CI_upper, fill = model ), alpha = 0.2 ) +
      #scale_color_brewer(palette = "Set1") +
      #scale_fill_brewer(palette = "Set1") +
      scale_color_manual(values = c("black", "blue")) +
      scale_fill_manual(values = c("black", "blue")) +
      geom_point( x = fd_exp_mean , y = fd_outcome_mean , size = 3, color = "darkorange1", shape = 15 ) +
      labs(title = "MR relationship", x = paste0("genotype predicted ", exposure),
           subtitle = paste0( "parametic-v-smooth: F-test P = ", MR_nonlinear_P, "; LRT P = ", MR_nonlinear_LRT_P) ) +
      theme_bw() +
      theme(legend.position = "none")
  }



  #############################
  ## IV. PLOT strata point
  ##     estimates
  #############################
  ###############################
  ## Stratified OBS. Estimates Plot
  ###############################
  if(plot_obs_res_betas == TRUE){ plotdata = obs_res_strata_ss} else { plotdata = obs_strata_ss }
  l_est_plot = plotdata %>%
    ggplot(aes(x = exposure_mean, y = beta)) +
    scale_fill_brewer(palette = brewer_col) +
    geom_errorbar( aes(ymin=beta-(1.96*se), ymax=beta+(1.96*se)), width = 0.2 ) +
    geom_point(aes(color = sig, shape = shape ) , size = 3) +
    scale_color_manual(values = c("black","red"), drop = FALSE) +
    scale_shape_manual(values = c(19,22) ) +
    geom_hline(yintercept = 0, size = 1, color = "black", linetype = "dashed" ) +
    theme_bw()

  if(plot_obs_res_betas == TRUE){
  l_est_plot = l_est_plot + labs(x = "mean of exposure strata", title = "GAM res. observational strata estimates",
       subtitle = paste0( "heterogenity P = ", res_obs_het_P, "; trend P = ", res_obs_trend_P) )
  } else {
    l_est_plot = l_est_plot + labs(x = "mean of exposure strata", title = "observational strata estimates",
                                   subtitle = paste0( "heterogenity P = ", res_obs_het_P, "; trend P = ", res_obs_trend_P) )
  }
  ## ADD strata color boundries
  if(add_strata_2_points == TRUE){
    l_est_plot = l_est_plot + geom_rect( aes(xmin = exposure_min, xmax = exposure_max,
                                            ymin = -Inf, ymax = Inf, fill = strataID), alpha = 0.2)
  }

  ###############################
  ## Stratified MR Estimates Plot
  ###############################
  iv_est_plot = tsls_strata_ivreg_ss %>%
    ggplot(aes(x = exposure_mean, y = beta )) +
    scale_fill_brewer(palette = brewer_col) +
    geom_errorbar( aes(ymin=beta-(1.96*se), ymax=beta+(1.96*se)), width = 0.2 ) +
    geom_point(aes(color = sig, shape = shape  ) , size = 3) +
    scale_color_manual(values = c("black","red"), drop = FALSE) +
    scale_shape_manual(values = c(19,22) ) +
    geom_hline(yintercept = 0, size = 1, color = "black", linetype = "dashed" ) +
    labs(x = "mean of exposure strata", y = "beta" ,title = "tsls MR strata estimates",
         subtitle = paste0( "heterogenity P = ", MR_het_P, "; trend P = ", MR_trend_P) ) +
    theme_bw()
  ## ADD strata color boundries
  if(add_strata_2_points == TRUE){
    iv_est_plot = iv_est_plot + geom_rect( aes(xmin = exposure_min, xmax = exposure_max,
                                               ymin = -Inf, ymax = Inf, fill = strataID), alpha = 0.2)
  }

  #############################
  ## V. MERGE PLOTS
  #############################
  if(old_plot_scheme == TRUE){
    plot1 = ggpubr::ggarrange(o_lnl_plot, iv_lnl_plot, nrow = 1, common.legend = TRUE,legend = "right" )
    plot2 = ggpubr::ggarrange(l_est_plot, iv_est_plot, nrow = 1, common.legend = TRUE,legend = "right" )
    plot = ggpubr::ggarrange(plot1, plot2, nrow = 2)
  } else {
    plot2 = ggpubr::ggarrange(l_est_plot, iv_est_plot, nrow = 1, common.legend = TRUE,legend = "right" )
    plot = ggpubr::ggarrange(o_lnl_plot, plot2, nrow = 1, widths = c(1,2))
  }

  #############################
  ## VI. Return to user
  #############################
  return(plot)

}
