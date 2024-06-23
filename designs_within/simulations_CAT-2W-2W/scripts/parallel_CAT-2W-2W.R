source("scripts/parallel_setup_script.R")

# CAT-2W-2W ---------------------------------------------------------------

cat_2w2w <- readr::read_csv("data/simulation_grids/CAT-2W-2W.csv")

cat_2w2w_combinations <- dplyr::full_join(
  tidyr::expand_grid(
    design = "CAT-2W-2W", 
    combo_id = cat_2w2w$combo_id,
    b, n_rm, 
    iter = start_iter:n_iter
  ), 
  cat_2w2w |> dplyr::select(-design),
  by = c("combo_id")
)  |> 
  dplyr::mutate(
    population_id = paste0(design, "_", combo_id, "_b_", b)
  )


generate_cat2w2w <- function(
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
    sd = 1, 
    r = 0.5, 
    within = list(w_group1 = c("gw1", "gw2"), w_group2 = c("a", "b")),
    plot = FALSE
  )
  
  transformed_df <- norm_df |> 
    dplyr::mutate(
      gw1_a = transform_norm(raw_var = gw1_a, mu = 0, nu = nu, tau = tau,
                             sigma = sigma),
      gw1_b = transform_norm(raw_var = gw1_b, mu = 0, nu = nu, tau = tau, 
                             sigma = sigma + sigma_shift_1) + b + b_shift_1, 
      gw2_a = transform_norm(raw_var = gw2_a, mu = 0, nu = nu, tau = tau,
                             sigma = sigma) + b,
      gw2_b = transform_norm(raw_var = gw2_b, mu = 0, nu = nu, tau = tau, 
                             sigma = sigma + sigma_shift_2) + 3*b + b_shift_2
    )
  
  return(transformed_df)
}

## sample fit and save ----------------------------------------------------

cat_2w2w_sfs <- function(
   # pop_i, 
    n_rm
    ){
  
  ## create a sample   
  sample_ij <- dplyr::slice_sample(.data = pop_i, n = n_rm)
  #sample_ij <- dplyr::slice_sample(.data = pop_i, n = 100)
  
  ## fit ols model 
  
  sample_ij_long <- sample_ij |> 
    tidyr::pivot_longer(-id, names_to = "gw", values_to = "y") |> 
    dplyr::mutate(
      gw2 = ifelse(stringr::str_detect(gw, "a"), "a", "b"), 
      gw1 = as.character(readr::parse_number(gw))
    )
  
  afex_summary_ij <- afex_nc_gg_summary("y ~ gw1*gw2 + (gw1*gw2|id)", data = sample_ij_long)
  
  rob_trim_summary_ij <- WRS::wwtrim(x = sample_ij[-1], J = 2, K = 2, tr = 0.2)
  rob_trim_summary_ij <- tibble::tibble(
    statistic = unlist(rob_trim_summary_ij[c(1,3,5)]), 
    p = unlist(rob_trim_summary_ij[c(2,4,6)])
  )
  names(rob_trim_summary_ij) <- paste0("rob_trim_", names(rob_trim_summary_ij))
  
  rob_boot_summary_ij <- WRS::wwtrimbt(x = sample_ij[-1], J = 2, K = 2, tr = 0) 
  rob_boot_summary_ij <- tibble::tibble(
    rob_boot_p = unlist(rob_boot_summary_ij)
  )
  
  rob_trimboot_summary_ij <- WRS::wwtrimbt(x = sample_ij[-1], J = 2, K = 2, tr = 0.2) 
  rob_trimboot_summary_ij <- tibble::tibble(
    rob_trimboot_p = unlist(rob_trimboot_summary_ij)
  )
  
  afex_int_b <- WRS::wwmcp(x = sample_ij[-1], J = 2, K = 2, tr = 0)$Factor_AB$psihat[2]
  rob_trim_int_b <- WRS::wwmcp(x = sample_ij[-1], J = 2, K = 2, tr = 0.2)$Factor_AB$psihat[2]
  
  dplyr::bind_cols(
    afex_summary_ij, rob_trim_summary_ij, rob_boot_summary_ij, rob_trimboot_summary_ij, 
    afex_int_b = afex_int_b, rob_trim_int_b = rob_trim_int_b
  )
  
  
  
}


## Iterate -----------------------------------------------------------------

set.seed(global_seed)

cat_2w2w_unique_populations <- cat_2w2w_combinations |> 
  dplyr::filter(!duplicated(population_id))

cat_2w2w_n_pop <- nrow(cat_2w2w_unique_populations)

for(i in 51:cat_2w2w_n_pop){
  
  show_progress(i = i, n_pop = cat_2w2w_n_pop)
  
  pop_row_i <- cat_2w2w_unique_populations[i, ]
  pop_i <- generate_cat2w2w(
    n = n_population, b = pop_row_i$b, 
    sigma_shift_1 = pop_row_i$sigma_shift_1, sigma_shift_2 = pop_row_i$sigma_shift_2, 
    nu = pop_row_i$nu, tau = pop_row_i$tau
  )
  
  grid_i <- cat_2w2w_combinations |> dplyr::filter(population_id == pop_row_i$population_id)
  
  ## start furrr level: sample fit and save
  
  tictoc::tic()
  future::plan(multisession)
  
  furrr_seed <- ifelse((start_iter == 1 & n_iter == 500), global_seed-1, global_seed)
  
  sim_results = furrr::future_pmap(
    .l = grid_i |> dplyr::select(n_rm), 
    .f = cat_2w2w_sfs,
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

