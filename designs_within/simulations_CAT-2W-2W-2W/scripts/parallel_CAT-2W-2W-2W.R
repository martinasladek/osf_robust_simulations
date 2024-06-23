source("scripts/parallel_setup_script.R")


# CAT-2W-2W-2W ------------------------------------------------------------

cat_2w2w2w <- readr::read_csv("data/simulation_grids/CAT-2W-2W-2W.csv")

cat_2w2w2w_combinations <- dplyr::full_join(
  tidyr::expand_grid(
    design = "CAT-2W-2W-2W", 
    combo_id = cat_2w2w2w$combo_id,
    b, n_rm, 
    iter = start_iter:n_iter
  ), 
  cat_2w2w2w |> dplyr::select(-design),
  by = c("combo_id")
) |> 
  dplyr::mutate(
    population_id = paste0(design, "_", combo_id, "_b_", b)
  )


generate_cat2w2w2w <- function(
    n = 100000, 
    b, 
    b_shift_1 = 0, 
    b_shift_2 = 0, 
    sigma = 1, 
    sigma_shift_1, 
    sigma_shift_2,
    nu, 
    tau
){
  
  set.seed(global_seed) 
  norm_df <- faux::sim_design(
    n = n, 
    mu = 0, 
    within = list(c("gw1", "gw2"), factor_2 = c("a", "b"), factor_3 = c(".", "..")),
    plot = FALSE
  )
  
  transformed_df <- norm_df |> 
    dplyr::mutate(
      `gw1_a_.` = transform_norm(raw_var = `gw1_a_.`, mu = 0, nu = nu, tau = tau, sigma = sigma),
      `gw1_a_..` = transform_norm(raw_var = `gw1_a_..`, mu = 0, nu = nu, tau = tau, sigma = sigma) + b,
      `gw1_b_.` = transform_norm(raw_var = `gw1_b_.`, mu = 0, nu = nu, tau = tau, sigma = sigma) + b,
      `gw1_b_..` = transform_norm(raw_var = `gw1_b_..`, mu = 0, nu = nu, tau = tau, sigma = sigma + sigma_shift_1) + 3*b + b_shift_1,
      
      `gw2_a_.` = transform_norm(raw_var = `gw2_a_.`, mu = 0, nu = nu, tau = tau, sigma = sigma) + b,
      `gw2_a_..` = transform_norm(raw_var = `gw2_a_..`, mu = 0, nu = nu, tau = tau, sigma = sigma) + 3*b,
      `gw2_b_.` = transform_norm(raw_var = `gw2_b_.`, mu = 0, nu = nu, tau = tau, sigma = sigma) + 3*b,
      `gw2_b_..` = transform_norm(raw_var = `gw2_b_..`, mu = 0, nu = nu, tau = tau, sigma = sigma + sigma_shift_2) + 7*b + b_shift_2,
    )
  
  return(transformed_df)
}


## sample fit and save ------------------------------------------------------
cat_2w2w2w_sfs <- function(
    #pop_i, 
    n_rm){
  
  ## create a sample   
  sample_ij <- dplyr::slice_sample(.data = pop_i, n = n_rm)
  #sample_ij <- dplyr::slice_sample(.data = pop_i, n = 18)
  
  ## fit ols model 
  
  sample_ij_long <- sample_ij |> 
    tidyr::pivot_longer(-id, names_to = "gw", values_to = "y") |> 
    dplyr::mutate(
      gw1 = as.character(readr::parse_number(gw)),
      gw2 = ifelse(stringr::str_detect(gw, "a"), "a", "b"), 
      gw3 = ifelse(stringr::str_detect(gw, "\\.{2}"), "..", ".")
    ) 
  
  afex_summary_ij <- afex_nc_gg_summary("y ~ gw1*gw2*gw3 + (gw1*gw2*gw3|id)", data = sample_ij_long)
  
  ## fit robust model(s)
  
  rob_trim_mod_ij <- WRS:::wwwtrim(2,2,2, sample_ij[-1], tr = 0.2)
  rob_trim_summary_ij <- tibble::tibble(
    rob_trim_statistic = unlist(rob_trim_mod_ij[grep(pattern = "p.value",  names(rob_trim_mod_ij), invert = TRUE)]), 
    rob_trim_p = unlist(rob_trim_mod_ij[grep(pattern = "p.value",  names(rob_trim_mod_ij), invert = FALSE)])
  )
  
  rob_boot_summary_ij <- tibble::tibble(
    boot_p = unlist(WRS::wwwtrimbt(2,2,2, sample_ij[-1], tr = 0, SEED = FALSE))
  )
  
  rob_trimboot_summary_ij <- tibble::tibble(
    trimboot_p = unlist(WRS::wwwtrimbt(2,2,2, sample_ij[-1], tr = 0.2, SEED = FALSE))
  )
  
  afex_int_b <- WRS::rm3mcp(2,2,2, sample_ij[-1], tr = 0)$Factor.ABC$psihat[2]
  rob_trim_int_b <- WRS::rm3mcp(2,2,2, sample_ij[-1], tr = 0.2)$Factor.ABC$psihat[2]
  
  dplyr::bind_cols(
    afex_summary_ij, rob_trim_summary_ij, rob_boot_summary_ij, rob_trimboot_summary_ij, 
    afex_int_b = afex_int_b, rob_trim_int_b = rob_trim_int_b
  )
  
}


## Iterate -----------------------------------------------------------------

set.seed(global_seed)

cat_2w2w2w_unique_populations <- cat_2w2w2w_combinations |> 
  dplyr::filter(!duplicated(population_id))

cat_2w2w2w_n_pop <- nrow(cat_2w2w2w_unique_populations)

cat_2w2w2w_sim_results <- data.frame()

for(i in 52:cat_2w2w2w_n_pop){
  
  show_progress(i = i, n_pop = cat_2w2w2w_n_pop)
  
  pop_row_i <- cat_2w2w2w_unique_populations[i, ]
  pop_i <- generate_cat2w2w2w(
    n = n_population, b = pop_row_i$b, 
    sigma_shift_1 = pop_row_i$sigma_shift_1, sigma_shift_2 = pop_row_i$sigma_shift_2, 
    nu = pop_row_i$nu, tau = pop_row_i$tau
  )
  
  grid_i <- cat_2w2w2w_combinations |> dplyr::filter(population_id == pop_row_i$population_id)
  
  ## start furrr level: sample fit and save
  
  tictoc::tic()
  future::plan(multisession)
  
  furrr_seed <- ifelse((start_iter == 1 & n_iter == 500), global_seed-1, global_seed)
  
  sim_results = furrr::future_pmap(
    .l = grid_i |> dplyr::select(n_rm), 
    .f = cat_2w2w2w_sfs,
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
