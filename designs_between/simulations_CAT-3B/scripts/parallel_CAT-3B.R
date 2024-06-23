
# CAT-3B  -----------------------------------------------------------------

source("scripts/parallel_setup_script.R")

cat_3b <- readr::read_csv("data/simulation_grids/CAT-3B.csv")

cat_3b_combinations <- dplyr::full_join(
  tidyr::expand_grid(
    design = "CAT-3B", 
    combo_id = cat_3b$combo_id,
    b, n_bw, 
    n_ratio = n_ratio_bw, 
    iter = 1:n_iter
  ), 
  cat_3b |> dplyr::select(-design),
  by = c("combo_id")
) |> 
  dplyr::mutate(
    population_id = paste0(design, "_", combo_id, "_b_", b)
  )


# Generate population -----------------------------------------------------

generate_cat3b <- function(
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
    between = list(gb1 = c("a", "b", "c")), 
    plot = FALSE
  )
  
  transformed_df <- norm_df |> 
    dplyr::mutate(
      y = dplyr::case_when(
        gb1 == "a" ~ transform_norm(raw_var = y, mu = 0, sigma = sigma, nu = nu, tau = tau),
        gb1 == "b" ~ transform_norm(raw_var = y, mu = 0, sigma = (sigma), nu = nu, tau = tau) + b,
        gb1 == "c" ~ transform_norm(raw_var = y, mu = 0, sigma = (sigma+sigma_shift), nu = nu, tau = tau) + (b + b_shift)
      )
    )
  
  return(transformed_df)
}



# sample fit and save -----------------------------------------------------


cat_3b_sfs <- function(n_bw, n_ratio){
  
  ## calculate sample split 
  ns <- sample_n_from_ratio_bw3(n_bw, n_ratio)
  
  ## create a sample   
  sample_ij <- dplyr::bind_rows(
    pop_i_a |> dplyr::slice_sample(n = ns[["n1"]]), 
    pop_i_b |> dplyr::slice_sample(n = ns[["n2"]]), 
    pop_i_c |> dplyr::slice_sample(n = ns[["n3"]]) 
  )
  
  # Fit models and export: 
  fit_and_summarise_mods_bw(mod_formula = "y~gb1", sample_ij = sample_ij)
  
}


# Iterate -----------------------------------------------------------------

set.seed(global_seed)

cat_3b_unique_populations <- cat_3b_combinations |> 
  dplyr::filter(!duplicated(population_id))

cat_3b_n_pop <- nrow(cat_3b_unique_populations)

for(i in 1:cat_3b_n_pop){
  
  show_progress(i = i, n_pop = cat_3b_n_pop)
  
  pop_row_i <- cat_3b_unique_populations[i, ]
  pop_i <- generate_cat3b(
    n = n_population, b = pop_row_i$b, 
    sigma_shift = pop_row_i$sigma_shift, 
    nu = pop_row_i$nu, tau = pop_row_i$tau
  )
  
  pop_i_a <- pop_i |> dplyr::filter(gb1 == "a") 
  pop_i_b <- pop_i |> dplyr::filter(gb1 == "b") 
  pop_i_c <- pop_i |> dplyr::filter(gb1 == "c") 
  
  grid_i <- cat_3b_combinations |> dplyr::filter(population_id == pop_row_i$population_id)
  
  
  ## start furrr level: sample fit and save
  
  tictoc::tic()
  future::plan(multisession)
  
  sim_results = furrr::future_pmap(
    .l = grid_i |> dplyr::select(n_bw, n_ratio), 
    .f = cat_3b_sfs,
    .options = furrr::furrr_options(seed = global_seed),
    .progress = TRUE
  )
  
  grid_i <- dplyr::bind_cols(
    grid_i,
    tibble::tibble(sim_results)
  )|> tidyr::unnest(sim_results)
  
  # Export
  export_results(grid_i, i=i, pop_row_i$population_id)
  tictoc::toc()
  
  rm(pop_i_a)
  rm(pop_i_b)
  rm(pop_i_c)
  
}


rm(pop_row_i)
rm(pop_i)
rm(grid_i)
gc()
