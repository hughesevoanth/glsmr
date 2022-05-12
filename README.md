# glsmr

### GAM and linear stratified MR

This is an R package to aid in determining if observational or two-stage least square (instrument variable or in genetics Mendelian randomization) analysis have a non-linear relationship between expoosure and outcome. 

### installation

	devtools::install_github("hughesevoanth/glsmr")

### use
	
- There are two functions that are most useful
	- glsmr()
	- plot_glsmr()
	
- an example for using glsmr
       	
		myexample = glsmr( wdata = mydata,
	          outcome = "trait",
	          exposure = "bmi",
	          instrument = "bmi_grs",
	          linear_covariates = c("batch", "sex"),
	          smooth_covariates =  c("age"),
	          # strata = 4, ## for quartiles
	          strata = c(10,18.5,25,30,45),
	          rnt_outcome = TRUE,
	          weights_variable = NA,
	          outlier_method = "iqr",
	          outlier_cutoff = 5,
	          messages = FALSE,
	          return_models = TRUE)
	

	- NOTE: 'smooth_covariates' will be modeled as smooths or non-linear variables in the GAM, but as typical parametric variables in the linear and tsls analyses. 
          
- an example for using plot_glsmr()
		
		plot_glsmr(myexample,
			add_strata_2_curves = FALSE,
			add_strata_2_points = TRUE,
			brewer_col = "Set1",
			old_plot_scheme = TRUE,
			old_GAM_smooths = TRUE,
			plot_obs_res_betas = FALSE,
		  	pval_thresh = 0.05)


### An example figure from plot_glsmr()

![](figures/plot_v3.png)

### See this repo's wiki for more information

[wiki page here](https://github.com/hughesevoanth/glsmr/wiki)



