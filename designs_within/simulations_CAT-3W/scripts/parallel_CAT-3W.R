source("scripts/parallel_setup_script.R")


# CAT-3W ------------------------------------------------------------------


cat_3w <- readr::read_csv("data/simulation_grids/CAT-3W.csv")

cat_3w_combinations <- dplyr::full_join(
  tidyr::expand_grid(
    design = "CAT-3w", 
    combo_id = cat_3w$combo_id,
    b, n_rm, 
    iter = start_iter:n_iter
  ), 
  cat_3w |> dplyr::select(-design),
  by = c("combo_id")
) |> 
  dplyr::mutate(
    population_id = paste0(design, "_", combo_id, "_b_", b)
  )


generate_cat3w <- function(
    n = 100000, 
    b, 
    b_shift = 0, 
    sigma = 1, 
    sigma_shift,
    nu, 
    tau
){
  
  set.seed(global_seed) 
  norm_df <- faux::sim_design(
    n = n, 
    mu = 0, 
    sd = 1, 
    r = 0.5, 
    within = list(gw1 = c("gw1", "gw2", "gw3")), 
    plot = FALSE
  )
  
  transformed_df <- norm_df |> 
    dplyr::mutate(
      gw1 = transform_norm(raw_var = gw1, mu = 0, sigma = sigma, nu = nu, tau = tau),
      gw2 = transform_norm(raw_var = gw2, mu = 0, sigma = sigma, nu = nu, tau = tau) + b,
      gw3 = transform_norm(raw_var = gw3, mu = 0, sigma = sigma + sigma_shift, nu = nu, tau = tau) + b
    )
  
  return(transformed_df)
}

## sample fit and save ------------------------------------------------------

cat_3w_sfs <- function(n_rm){
  
  ## create a sample   
  sample_ij <- dplyr::slice_sample(.data = pop_i, n = n_rm)
  #sample_ij <- dplyr::slice_sample(.data = pop_i, n = 100)
  
  sample_ij_long <- sample_ij |> 
    tidyr::pivot_longer(-id, names_to = "gw", values_to = "y") 
  
  ## fit ols model 
  
  afex_nc_gg_summary_ij <- afex_nc_gg_summary("y ~ gw + (gw|id)", data = sample_ij_long) 
  afex_params_ij <- dplyr::bind_rows(
    afex_ttest_params("gw2", "gw1", sample_ij, term = "gw2-gw1"), 
    afex_ttest_params("gw3", "gw1", sample_ij, term = "gw3-gw1")
  )
  
  rob_trim_mod_ij <- WRS::rmanova(sample_ij[, -1])
  rob_trim_summary_ij <- format_wrs_mod(rob_trim_mod_ij, c("test", "p.value"), "rob_trim")
  rob_trim_params_ij <- dplyr::bind_rows(
    rob_trim_params("gw2", "gw1", sample_ij, term = "gw2-gw1"), 
    rob_trim_params("gw3", "gw1", sample_ij, term = "gw3-gw1")
  )
  
  rob_boot_mod_ij <- WRS::rmanovab(x = sample_ij[, -1], tr = 0) |> stfu()
  rob_boot_summary_ij <- format_wrs_mod(rob_boot_mod_ij, c("teststat", "p.value"), "rob_boot")
  rob_boot_params_ij <- dplyr::bind_rows(
    rob_trimboot_params("gw2", "gw1", sample_ij, term = "gw2-gw1", tr = 0, prefix = "boot"), 
    rob_trimboot_params("gw3", "gw1", sample_ij, term = "gw3-gw1", tr = 0, prefix = "boot")
  )
  
  rob_trimboot_mod_ij <- WRS::rmanovab(x = sample_ij[, -1], tr = 0.2) |> stfu()
  rob_trimboot_summary_ij <- format_wrs_mod(rob_trimboot_mod_ij, c("teststat", "p.value"), "rob_trimboot")
  rob_trimboot_params_ij <- dplyr::bind_rows(
    rob_trimboot_params("gw2", "gw1", sample_ij, term = "gw2-gw1", tr = 0.2, prefix = "trimboot"), 
    rob_trimboot_params("gw3", "gw1", sample_ij, term = "gw3-gw1", tr = 0.2, prefix = "trimboot")
  )
  
  dplyr::bind_cols(
    afex_nc_gg_summary_ij, afex_params_ij, 
    rob_trim_summary_ij, rob_trim_params_ij, 
    rob_boot_summary_ij, rob_boot_params_ij, 
    rob_trimboot_summary_ij, rob_trimboot_params_ij
  )
  
}

## Iterate -----------------------------------------------------------------

  
set.seed(global_seed)

cat_3w_unique_populations <- cat_3w_combinations |> 
  dplyr::filter(!duplicated(population_id))

cat_3w_n_pop <- nrow(cat_3w_unique_populations)

for(i in 34:cat_3w_n_pop){
  
  show_progress(i = i, n_pop = cat_3w_n_pop)
  
  pop_row_i <- cat_3w_unique_populations[i, ]
  pop_i <- generate_cat3w(
    n = n_population, b = pop_row_i$b, 
    sigma_shift = pop_row_i$sigma_shift, 
    nu = pop_row_i$nu, tau = pop_row_i$tau
  )
  
  grid_i <- cat_3w_combinations |> dplyr::filter(population_id == pop_row_i$population_id)
  
  
  ## start furrr level: sample fit and save
  
  tictoc::tic()
  future::plan(multisession)
  
  furrr_seed <- ifelse((start_iter == 1 & n_iter == 500), global_seed-1, global_seed)
  
  sim_results = furrr::future_pmap(
    .l = grid_i |> dplyr::select(n_rm), 
    .f = cat_3w_sfs,
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



