gc()

# Packages --------------------------------------------------------------------
library(afex)
library(dplyr)
library(faux)
library(furrr)
library(ggplot2)
library(parameters)
library(purrr)
library(tictoc)
library(WRS)
library(WRS2)
library(zoo)

here::here("scripts/helpers.R") |> source()

# Parameter settings ----------------------------------------------------------

check_iterations <- function() {
  repeat {
    user_input <- readline("Have you set the correct number of iterations? (y/n): ")
    # Check the user input and act accordingly
    if (tolower(user_input) == "y") {
      print("Continuing...")
      break  # Exit the loop and continue with the script
    } else if (tolower(user_input) == "n") {
      print("Exiting...")
      message("Press 'y' to continue.")
    } else {
      print("Invalid input. Please enter 'y' or 'n'.")
    }
  }
}

#check_iterations()

start_iter <- 1
n_iter = 1000
n_population = 100000
global_seed <- 21023

# Shared population parameters 
b <- c(0, 0.1, 0.3, 0.5)

n_ratio_bw <- c(1, 1.5, 2)
n_ratio_rm <- c(1, 1.2, 2)

# Repeated measures sample metrics 
n_rm <- c(18, 65, 172)

# Between measures sample metrics
n_bw <- c(71, 296, 1141)

# Functions -------------------------------------------------------------------

check_iterations <- function() {
  repeat {
    user_input <- readline("Have you set the correct number of iterations? (y/n): ")
    # Check the user input and act accordingly
    if (tolower(user_input) == "y") {
      print("Continuing...")
      break  # Exit the loop and continue with the script
    } else if (tolower(user_input) == "n") {
      print("Exiting...")
      message("Press 'y' to continue.")
    } else {
      print("Invalid input. Please enter 'y' or 'n'.")
    }
  }
}

show_progress <- function(i, n_pop){
  cat(
    paste0("Population ", i, " out of ", n_pop, " (", i/n_pop*100, "%).\n\n")
  )
}

export_results <- function(result_file, i,  design_label){
  
  file_name <- paste0(
    design_label,
    "_i", i, 
    "_", Sys.Date(), ".rds"
  )
  
  file_path <- paste0("data/simulation_exports/", file_name)
  
  saveRDS(result_file, file_path)
  
}

transform_norm <- function(raw_var, mu, sigma, nu, tau){
  
  pnorm_var <- pnorm(raw_var)
  transformed_var <- gamlss.dist::qSHASHo(
    p = pnorm_var, mu = mu, sigma = sigma, nu = nu, tau = tau)
  
  return(transformed_var)
  
}

sample_n_from_ratio_rm <- function(n_rm, n_ratio){
  
  n1 = floor(n_rm / (n_ratio + 1))
  n2 = n_rm - n1
  
  return(
    c(n1 = n1, n2 = n2)
  )
}

sample_n_from_ratio_bw2 <- function(n_bw, n_ratio){
  
  n1 = floor(n_bw / (n_ratio + 1))
  n2 = n_bw - n1
  
  return(
    c(n1 = n1, n2 = n2)
  )
}

sample_n_from_ratio_bw3 <- function(n_bw, n_ratio) {
  
  n1 = n_bw / (2 + n_ratio)
  n2 = n1
  n3 = n_ratio * n1
  
  return(
    floor(c(n1 = n1, n2 = n2, n3 = n3))
  )
}

vp1 <- function(x1, vp_val) 1
vp23 <- function(x1, vp_val) (5 / (1+exp(-vp_val*(abs(x1)-1))))/2.5
vp4 <- function(x1, vp_val) (5 / (1+exp(-vp_val*((x1)-1))))/2.5

squish <- function(x, new_min = min_bound, new_max = max_bound){   
  (x - min(x))/(max(x)-min(x)) * (new_max - new_min) + new_min 
}

mod_fun_vp1 <- function(x) (x)

mod_fun_vp23 <- function(x){
  abs(x)
}

mod_fun_vp4 <- function(x){
  min = 0
  max = 1.75
  (max-min) * ((x - min(x)) / (max(x) - min(x)) + min)
}

text_as_fun <- function(function_string){
  eval(parse(text = function_string))
}

lm_summary <- function(lm_mod_ij){
  
  lm_mod_ij_betas <- broom::tidy(lm_mod_ij, conf.int = TRUE)
  names(lm_mod_ij_betas) <- paste0("lm_", names(lm_mod_ij_betas))
  
  lm_mod_ij_overall <- 
    broom::glance(lm_mod_ij, conf.int = TRUE)[, c("r.squared", "adj.r.squared", "statistic", "p.value")]
  names(lm_mod_ij_overall) <- paste0("lm_overall_", names(lm_mod_ij_overall))
  
  dplyr::bind_cols(lm_mod_ij_betas, lm_mod_ij_overall)
  
}

rob_mm_summary <- function(rob_mm_mod_ij, prefix){
  
  rob_mm_confint_ij <- tibble::as_tibble(confint(rob_mm_mod_ij))
  names(rob_mm_confint_ij) <- c("conf.low", "conf.high")
  rob_mm_ij_overall <- summary(rob_mm_mod_ij)
  
  rob_mm_mod_ij_betas <- dplyr::bind_cols(
    broom::tidy(rob_mm_mod_ij), rob_mm_confint_ij, 
    r.squared = rob_mm_ij_overall$r.squared, 
    adj.r.squared = rob_mm_ij_overall$adj.r.squared
  ) |> 
    dplyr::mutate(
      warning = NA
    )
  
  names(rob_mm_mod_ij_betas) <- paste0("rob_", prefix, "_", names(rob_mm_mod_ij_betas))
  
  rob_mm_mod_ij_betas
  
}

rob_hc4_summary <- function(lm_mod_ij){
  tibble::tibble(parameters::parameters(lm_mod_ij, vcov = "HC4")) |> 
    dplyr::select(
      rob_hc4_std.error = SE, 
      rob_hc4_conf.low = CI_low, rob_hc4_conf.high = CI_high, 
      rob_hc4_p.value = p
    )
}

rob_boot_summary <- function(lm_mod_ij){
  tibble::tibble(
    parameters::parameters(lm_mod_ij, bootstrap = TRUE, ci_method = "bcai", iteractions = 599)
  ) |> 
    dplyr::select(
      rob_boot_estimate = Coefficient, 
      rob_boot_conf.low = CI_low, rob_boot_conf.high = CI_high, 
      rob_boot_p.value = p
    )
}

blank_robust_export <- function(){
  
  tibble::tibble(
    term = NA,
    estimate = NA , std.error = NA, statistic = NA, p.value = NA, conf.low = NA, conf.high = NA,
    r.squared = NA, adj.r.squared = NA
  )
  
}

fit_and_summarise_mods_bw <- function(mod_formula, sample_ij, mm_label = "mm", ks_label = "ks"){
  
  #mod_formula = lm_formula
  mod_formula <- as.formula(mod_formula)
  
  lm_mod_ij <- lm(mod_formula, data = sample_ij)
  lm_mod_ij$call$formula <- mod_formula #explicitly define formula, otherwise parameters::parameters cry. 
  
  lm_summary_ij <- lm_summary(lm_mod_ij = lm_mod_ij)
  
  
  ## fit robust model(s)
  
  ### MM estimator
  
  rob_mm_mod_ij <- tryCatch(
    robustbase::lmrob(mod_formula, data = sample_ij),
    warning = function(w) return(as.character(w))
  )
  
  if(is.character(rob_mm_mod_ij)){
    
    rob_mm_summary_ij <- blank_robust_export() |> 
      dplyr::mutate(warning = rob_mm_mod_ij)
    
    names(rob_mm_summary_ij) <- paste0("rob_", mm_label, "_", names(rob_mm_summary_ij))
    
  } else {
    
    rob_mm_summary_ij <- rob_mm_summary(rob_mm_mod_ij = rob_mm_mod_ij, mm_label)
    
  }
  
  ### KS2014 estimator
  
  rob_ks_mod_ij <- tryCatch(
    robustbase::lmrob(mod_formula, data = sample_ij, setting = "KS2014"),
    warning = function(w) return(as.character(w))
  )
  
  if(is.character(rob_ks_mod_ij)){
    
    rob_ks_summary_ij <- blank_robust_export() |> 
      dplyr::mutate(warning = rob_ks_mod_ij)
    
    names(rob_ks_summary_ij) <- paste0("rob_", ks_label, "_", names(rob_ks_summary_ij))
    
  } else {
    
    rob_ks_summary_ij <- rob_mm_summary(rob_mm_mod_ij = rob_ks_mod_ij, ks_label)
    
  }
  
  ### HC4 
  rob_hc4_summary_ij <- rob_hc4_summary(lm_mod_ij = lm_mod_ij)
  
  ### Bootstrap
  rob_boot_summary_ij <- rob_boot_summary(lm_mod_ij = lm_mod_ij)
  
  
  # Return: 
  tibble::tibble(
    lm_summary = list(lm_summary_ij), 
    rob_mm_summary = list(rob_mm_summary_ij), 
    rob_ks_summary = list(rob_ks_summary_ij), 
    rob_hc4_summary = list(rob_hc4_summary_ij), 
    rob_boot_summary = list(rob_boot_summary_ij)
  ) |> 
    tidyr::unnest(cols = everything())
  
}

afex_summary <- function(afex_mod, prefix){
  
  afex_summary <- as.data.frame(parameters::parameters(afex_mod)) |> 
    dplyr::select(term = Parameter, df1 = df, df2 = df_error, statistic = `F`, p)
  
  names(afex_summary) <- paste0("afex_", prefix, "_", names(afex_summary))
  
  afex_summary
  
}

afex_nc_gg_summary <- function(afex_formula, data, factorize = TRUE){
  
  afex_formula = as.formula(afex_formula)
  
  suppressWarnings({
    
    afex_nc_mod <- afex::aov_4(afex_formula, data = data, factorize = factorize, anova_table = list(correction = "none")) 
    afex_nc_summary <- afex_summary(afex_nc_mod, "nc")
    
    afex_gg_mod <- afex::aov_4(afex_formula, data = data, factorize = factorize, anova_table = list(correction = "GG"))
    afex_gg_summary <- afex_summary(afex_nc_mod, "gg")
    
  })
  
  dplyr::bind_cols(afex_nc_summary, afex_gg_summary)
  
}

afex_ttest_params <- function(group2, group1, df, term = "gw"){
  
  ttest_mod_ij <- t.test(df[[group2]], df[[group1]], paired = TRUE)
  ttest_summary_ij <- parameters::parameters(ttest_mod_ij) |> 
    tibble::as_tibble() |> 
    dplyr::transmute(term = term, diff = Difference, statistic = t, p, conf.low = CI_low, conf.high = CI_high)
  names(ttest_summary_ij) <- paste0("afex_", "t_", names(ttest_summary_ij))
  
  ttest_summary_ij
  
}

rob_trim_params <- function(group2, group1, df, term = "gw"){
  
  rob_trim_mod <- WRS::yuend(df[[group2]], df[[group1]])
  
  rob_trim_summary <- tibble::tibble(
    term = term, diff = rob_trim_mod$dif, statistic = rob_trim_mod$teststat, p = rob_trim_mod$p.value,
    conf.low = rob_trim_mod$ci[1], conf.high = rob_trim_mod$ci[2]
  )
  names(rob_trim_summary) <- paste0("rob_trim_t_", names(rob_trim_summary))
  
  rob_trim_summary
  
}

rob_trimboot_params <- function(group2, group1, df, tr, prefix, term = "gw"){
  
  rob_trimboot_mod <- WRS::ydbt(df[[group2]], df[[group1]], tr = tr, SEED = FALSE) 
  rob_trimboot_summary <- tibble::tibble(
    term = term, diff = rob_trimboot_mod$dif, p = rob_trimboot_mod$p.value,
    conf.low = rob_trimboot_mod$ci[1], conf.high = rob_trimboot_mod$ci[2]
  )
  names(rob_trimboot_summary) <- paste0("rob_", prefix, "_t_", names(rob_trimboot_summary))
  
  rob_trimboot_summary
  
}

format_wrs_mod <- function(wrs_mod, metrics, prefix){
  
  wrs_mod_summary <- as.data.frame(wrs_mod[metrics])
  names(wrs_mod_summary) <- paste0(prefix, "_", names(wrs_mod_summary))
  
  wrs_mod_summary
  
}

stfu <- function(...){
  SimDesign::quiet(...)
}

