
# CAT-2B-3B ---------------------------------------------------------------

start_pop_i <- 1

source("scripts/parallel_setup_script.R")

cat_2b3b <- readr::read_csv("data/simulation_grids/CAT-2B-3B.csv")

cat_2b3b_combinations <- dplyr::full_join(
  tidyr::expand_grid(
    design = "CAT-2B-3B", 
    combo_id = cat_2b3b$combo_id,
    b, n_bw, 
    n_ratio = n_ratio_bw, 
    iter = 1:n_iter
  ), 
  cat_2b3b |> dplyr::select(-design),
  by = c("combo_id")
) |> 
  dplyr::mutate(
    population_id = paste0(design, "_", combo_id, "_b_", b)
  )

generate_cat2b3b <- function(
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
    between = list(
      gb1 = c("1", "2"),
      gb2 = c("a", "b", "c")
    ),
    plot = FALSE
  )
  
  transformed_df <- norm_df |> 
    dplyr::mutate(
      y = dplyr::case_when(
        gb1 == "1" & gb2 == "a" ~ transform_norm(raw_var = y, mu = 0, nu = nu, tau = tau,
                                                 sigma = sigma),
        gb1 == "1" & gb2 == "b" ~ transform_norm(raw_var = y, mu = 0, nu = nu, tau = tau,
                                                 sigma = sigma) + b,
        gb1 == "1" & gb2 == "c" ~ transform_norm(raw_var = y, mu = 0, nu = nu, tau = tau,
                                                 sigma = sigma + sigma_shift_1) + b + b_shift_1,
        gb1 == "2" & gb2 == "a" ~ transform_norm(raw_var = y, mu = 0, nu = nu, tau = tau,
                                                 sigma = sigma) + b,
        gb1 == "2" & gb2 == "b" ~ transform_norm(raw_var = y, mu = 0, nu = nu, tau = tau,
                                                 sigma = sigma) + 2*b,
        gb1 == "2" & gb2 == "c" ~ transform_norm(raw_var = y, mu = 0, nu = nu, tau = tau,
                                                 sigma = sigma + sigma_shift_2) + 2*b + b_shift_2,
      )
    )
  
  return(transformed_df)
}



# sample fit and save -----------------------------------------------------


cat_2b3b_sfs <- function(n_bw, n_ratio){
  
  ns <- sample_n_from_ratio_bw2(n_bw = n_bw, n_ratio = n_ratio)
  
  weight_col <- paste0("n_ratio_w_", n_ratio)
  
  sample_ij <- dplyr::bind_rows(
    pop_i_gb1_1 |> dplyr::slice_sample(n = ns[["n1"]], weight_by = (!!sym(weight_col))), 
    pop_i_gb1_2 |> dplyr::slice_sample(n = ns[["n2"]], weight_by = (!!sym(weight_col))) 
  )
  
  n_cell_count <- sample_ij |> dplyr::count(gb1, gb2) |> dplyr::pull(n)
  
  # must sample at least 2 for each cell and all groups must be represented
  # in the cell grid 
  while (any(n_cell_count < 2) | length(n_cell_count) < 6) {
    
    sample_ij <- dplyr::bind_rows(
      pop_i_gb1_1 |> dplyr::slice_sample(n = ns[["n1"]], weight_by = (!!sym(weight_col))), 
      pop_i_gb1_2 |> dplyr::slice_sample(n = ns[["n2"]], weight_by = (!!sym(weight_col))) 
    )
    
    n_cell_count <- sample_ij |> dplyr::count(gb1, gb2) |> dplyr::pull(n)
    
  }
  
  sample_ij <- sample_ij |> dplyr::select(!dplyr::contains("n_ratio")) 
  
  # Fit models and export: 
  fit_and_summarise_mods_bw(mod_formula = "y ~ gb1+gb2", sample_ij = sample_ij)
  
}


# Iterate -----------------------------------------------------------------


set.seed(global_seed)

cat_2b3b_unique_populations <- cat_2b3b_combinations |> 
  dplyr::filter(!duplicated(population_id))

cat_2b3b_n_pop <- nrow(cat_2b3b_unique_populations)


for(i in start_pop_i:cat_2b3b_n_pop){
  
  show_progress(i = i, n_pop = cat_2b3b_n_pop)
  
  pop_row_i <- cat_2b3b_unique_populations[i, ]
  
  grid_i <- cat_2b3b_combinations |> dplyr::filter(population_id == pop_row_i$population_id)
  
  n_ratios <- sort(unique(grid_i$n_ratio))
  
  pop_i <- generate_cat2b3b(
    n = n_population, b = pop_row_i$b, 
    sigma_shift_1 = pop_row_i$sigma_shift_1, sigma_shift_2 = pop_row_i$sigma_shift_2, 
    nu = pop_row_i$nu, tau = pop_row_i$tau
  )
  
  pop_i[[paste0("n_ratio_w_", n_ratios[1])]] <- ifelse(pop_i$gb2 == "c" , n_ratios[1], 1)
  pop_i[[paste0("n_ratio_w_", n_ratios[2])]] <- ifelse(pop_i$gb2 == "c" , n_ratios[2], 1)
  pop_i[[paste0("n_ratio_w_", n_ratios[3])]] <- ifelse(pop_i$gb2 == "c" , n_ratios[3], 1)
  
  pop_i_gb1_1 <- pop_i |> dplyr::filter(gb1 == "1")
  pop_i_gb1_2 <- pop_i |> dplyr::filter(gb1 == "2")
  
  ## start furrr level: sample fit and save
  
  tictoc::tic()
  future::plan(multisession)
  
  sim_results = furrr::future_pmap(
    .l = grid_i |> dplyr::select(n_bw, n_ratio), 
    .f = cat_2b3b_sfs,
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
  
  rm(pop_i_gb1_1)
  rm(pop_i_gb1_2)
  
}

rm(pop_row_i)
rm(pop_i)
rm(grid_i)
gc()

```


