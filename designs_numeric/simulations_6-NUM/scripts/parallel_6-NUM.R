source("scripts/parallel_setup_script.R")

# N-NUM -------------------------------------------------------------------

n6_num <- read.csv("data/simulation_grids/6-NUM.csv")
n3_num <- read.csv("data/simulation_grids/3-NUM.csv")
n2_num <- read.csv("data/simulation_grids/2-NUM.csv")
n1_num <- read.csv("data/simulation_grids/1-NUM.csv")

n_num <- rbind(
  n6_num, n3_num, n2_num, n1_num
)

n_num_combinations <- dplyr::full_join(
  tidyr::expand_grid(
    design = "N-NUM",
    combo_id = n_num$combo_id,
    b, n_bw, n_preds = c(1,2,3,6),
    iter = start_iter:n_iter
  ), 
  n_num |> dplyr::select(-design), 
  by = c("combo_id", "n_preds")
) |> 
  dplyr::mutate(
    population_id = paste0(design, "_", combo_id, "_b_", b, "_n_preds_", n_preds)
  ) |> 
  dplyr::filter(
    !n_het_pred > n_preds
  )


generate_n_num <- function(
    n = 100000,
    b,
    n_preds,
    min_bound, 
    max_bound, 
    nu, 
    tau, 
    vp_val, 
    frac, 
    n_het_pred, 
    mod_fun,
    vp_fun
){
  
  varnames <- paste0("x", 1:n_preds)
  
  set.seed(global_seed) 
  norm_df <- faux::rnorm_multi(
    r = 0.212, 
    n = n, 
    mu = 0, 
    sd = 1, 
    varnames = varnames
  )
  
  # generate errors: 
  set.seed(global_seed) 
  e <- gamlss.dist::rSHASHo(
    n = n,
    mu = 0,
    nu = nu,
    tau = tau
  )
  
  if(n_het_pred == 1){
    
    x1_mod <- mod_fun(norm_df$x1)
    x_mult <- (x1_mod)^(1/frac)
    
    m_total <- vp_fun(x = x_mult, vp_val)
    
  } else if(n_het_pred == 2){
    
    x1_mod <- mod_fun(norm_df$x1)
    x2_mod <- mod_fun(norm_df$x2)
    x_mult <- (x1_mod*x2_mod)^(1/frac)
    
    m_total <- vp_fun(x = x_mult, vp_val)
    
  } else if(n_het_pred == 0) {
    
    m_total = 1
    
  }
  
  # modify the error: 
  e_mod <- (e*m_total)
  
  # adjust error to fit in between bounds 
  e_mod <- squish(e_mod, min_bound, max_bound)
  
  # compute y: 
  norm_df$y <- rowSums(norm_df*b, na.rm = TRUE) + e_mod
  
  return(norm_df)
}

## sample fit and save -----------------------------------------------------

n_num_sfs <- function(
    #pop_i, 
    n_bw){
  
  ## create a sample   
  sample_ij <- dplyr::slice_sample(.data = pop_i, n = n_bw)
  
  var_names <- paste0(names(sample_ij))
  var_names <- var_names[-length(var_names)]
  ## fit ols model 
  lm_formula <- paste0("y ~ ", paste0(var_names, collapse = " + "))
  
  # Fit models and export: 
  fit_and_summarise_mods_bw(mod_formula = lm_formula, sample_ij = sample_ij)
  
}


## Iterate -----------------------------------------------------------------

set.seed(global_seed)

n_num_unique_populations <- n_num_combinations |> 
  dplyr::filter(
    !duplicated(population_id), 
    n_preds %in% c(6)
  )

n_num_n_pop <- nrow(n_num_unique_populations)

for(i in 132:n_num_n_pop){
  
  show_progress(i = i, n_pop = n_num_n_pop)
  
  pop_row_i <- n_num_unique_populations[i, ]
  pop_i <- generate_n_num(
    n = n_population, b = pop_row_i$b, 
    n_preds = pop_row_i$n_preds, 
    min_bound = pop_row_i$min_bound, max_bound = pop_row_i$max_bound, 
    vp_val = pop_row_i$vp_val, frac = pop_row_i$frac, n_het_pred = pop_row_i$n_het_pred, 
    mod_fun = text_as_fun(pop_row_i$mod_fun), 
    vp_fun = text_as_fun(pop_row_i$vp_fun), 
    nu = pop_row_i$nu, tau = pop_row_i$tau
  )
  
  grid_i <- n_num_combinations |> dplyr::filter(population_id == pop_row_i$population_id)
  
  ## start furrr level: sample fit and save
  
  tictoc::tic()
  future::plan(multisession)
  
  furrr_seed <- ifelse((start_iter == 1 & n_iter == 500), global_seed-1, global_seed)
  
  sim_results = furrr::future_pmap(
    .l = grid_i |> dplyr::select(n_bw), 
    .f = n_num_sfs,
    .options = furrr::furrr_options(seed = furrr_seed),
    .progress = TRUE
  )
  
  grid_i <- dplyr::bind_cols(
    grid_i,
    tibble::tibble(sim_results)
  )|> tidyr::unnest(sim_results)
  
  # Export
  export_results(grid_i, i = paste0(i, "_", start_iter, "-", n_iter), pop_row_i$population_id)
  tictoc::toc()
  
  rm(pop_i)
  
}

rm(pop_row_i)
rm(pop_i)
rm(grid_i)
gc()
