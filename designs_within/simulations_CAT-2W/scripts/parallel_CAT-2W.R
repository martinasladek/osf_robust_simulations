source("scripts/parallel_setup_script.R")

# CAT-2W ------------------------------------------------------------------

cat_2w <- readr::read_csv("data/simulation_grids/CAT-2W.csv")

cat_2w_combinations <- dplyr::full_join(
  tidyr::expand_grid(
    design = "CAT-2W", 
    combo_id = cat_2w$combo_id,
    b, n_rm, 
    iter = start_iter:n_iter
  ), 
  cat_2w |> dplyr::select(-design),
  by = c("combo_id")
) |> 
  dplyr::mutate(
    population_id = paste0(design, "_", combo_id, "_b_", b)
  )


generate_cat2w <- function(
    n = 100000, 
    b, 
    b_shift = 0, 
    sigma_shift, 
    sigma = 1,
    nu, 
    tau
){
  
  set.seed(global_seed) 
  norm_df <- faux::sim_design(
    n = n, 
    mu = 0, 
    sd = 1, 
    r = 0.5, 
    within = list(w_group = paste0("gw", 1:2)), 
    plot = FALSE
  )
  
  transformed_df <- norm_df |> 
    dplyr::mutate(
      gw1 = transform_norm(raw_var = gw1, mu = 0, nu = nu, tau = tau,
                           sigma = sigma),
      gw2 = transform_norm(raw_var = gw2, mu = 0, nu = nu, tau = tau,
                           sigma = (sigma+sigma_shift)) + (b + b_shift)
    )
  
  return(transformed_df)
}


## sample fit and save -----------------------------------------------------

cat_2w_sfs <- function(
   # pop_i, 
    n_rm
){
  
  ## create a sample   
  sample_ij <- dplyr::slice_sample(.data = pop_i, n = n_rm)
  #sample_ij <- dplyr::slice_sample(.data = pop_i, n = 100)
  
  ## fit ols model 
  
  afex_summary_ij <- afex_ttest_params("gw2", "gw1", sample_ij)
  rob_trim_summary_ij <- rob_trim_params("gw2", "gw1", sample_ij)
  rob_boot_summary_ij <- rob_trimboot_params("gw2", "gw1", sample_ij, 0, "boot")
  rob_trimboot_summary_ij <- rob_trimboot_params("gw2", "gw1", sample_ij, 0.2, "trimboot")
  
  # Return: 
  dplyr::bind_cols(afex_summary_ij, rob_trim_summary_ij, rob_boot_summary_ij, rob_trimboot_summary_ij)
  
}


## Iterate -----------------------------------------------------------------

  
set.seed(global_seed)

cat_2w_unique_populations <- cat_2w_combinations |> 
  dplyr::filter(!duplicated(population_id))

cat_2w_n_pop <- nrow(cat_2w_unique_populations)

for(i in 2:cat_2w_n_pop){
  
  show_progress(i = i, n_pop = cat_2w_n_pop)
  
  pop_row_i <- cat_2w_unique_populations[i, ]
  pop_i <- generate_cat2w(
    n = n_population, b = pop_row_i$b, 
    sigma_shift = pop_row_i$sigma_shift, 
    nu = pop_row_i$nu, tau = pop_row_i$tau
  )
  
  grid_i <- cat_2w_combinations |> dplyr::filter(population_id == pop_row_i$population_id)
  
  ## start furrr level: sample fit and save
  
  tictoc::tic()
  future::plan(multisession)
  
  furrr_seed <- ifelse((start_iter == 1 & n_iter == 500), global_seed-1, global_seed)
  
  sim_results = furrr::future_pmap(
    .l = grid_i |> dplyr::select(n_rm), 
    .f = cat_2w_sfs,
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
