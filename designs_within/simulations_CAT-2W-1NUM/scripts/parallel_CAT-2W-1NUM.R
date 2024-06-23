source("scripts/parallel_setup_script.R")

# CAT-2W-1NUM -------------------------------------------------------------

cat_2w1num <- readr::read_csv("data/simulation_grids/CAT-2W-1NUM.csv")

cat_2w1num_combinations <- dplyr::full_join(
  tidyr::expand_grid(
    design = "CAT-2w1num", 
    combo_id = cat_2w1num$combo_id,
    b, n_rm, 
    iter = start_iter:n_iter
  ), 
  cat_2w1num |> dplyr::select(-design),
  by = c("combo_id")
) |> 
  dplyr::mutate(
    population_id = paste0(design, "_", combo_id, "_b_", b)
  )


generate_cat2w1num <- function(
    n, 
    b, 
    sigma = 1, 
    sigma_shift, 
    min_bound, 
    max_bound, 
    nu,
    tau, 
    vp_val, 
    vp_fun, 
    mod_fun
){
  
  #generate covariate: 
  cov = rnorm(n, 0, 1)
  
  #generate RM data: 
  set.seed(global_seed) 
  df <- faux::rnorm_multi(
    n = n, vars = 2, r = 0.5, varnames = c(0, 1)
  ) |> 
    dplyr::mutate(id = 1:(n)) |> 
    tidyr::pivot_longer(cols = -id, names_to = "gw",  values_to = "e") |> 
    dplyr::arrange(gw) |> 
    dplyr::mutate(
      gw = as.numeric(gw), 
      cov = c(cov, cov), 
      e = transform_norm(raw_var = e, mu = 0, sigma = sigma, nu = nu, tau = tau),
    )
  
  # calculate multipliers for setting up heteroscedasticty: 
  cov_mult <- mod_fun(df$cov)
  m_cont <- vp_fun(x = cov_mult, vp_val)
  m_cat = ifelse(df$gw == 1, sigma + sigma_shift, sigma)
  m_total = m_cont*m_cat
  #m_total = 1
  
  # adjust errors and compute y
  cov_shift = 0.25
  df <- df |> 
    dplyr::mutate(
      e = df$e * m_total,
      e = squish(e, min_bound, max_bound), 
      y = b*gw + (b-cov_shift)*cov + b*gw*cov + e
    ) |> 
    dplyr::select(-e)
  
  return(df)
  
}

## sample fit and save -------------------------------------------------------

cat_2w1num_sfs <- function(
    #pop_i, 
    n_rm){
  
  ## create a sample   
  sample_ij_wide <- pop_i_wide |> dplyr::slice_sample(n = n_rm)
  #sample_ij_wide <- pop_i_wide |> dplyr::slice_sample(n = 500)
  
  sample_ij <- sample_ij_wide |> tidyr::pivot_longer(cols = -c(id, cov), names_to = "gw", values_to = "y")
  
  #export_results(sample_ij, i = i, design_label = paste0("debugging"))

  
  afex_mod_ij <- afex::aov_4(y ~ gw*cov + (gw|id), data = sample_ij, factorize = FALSE) |> suppressWarnings()
  afex_summary_ij <- afex_nc_gg_summary("y ~ gw*cov + (gw|id)", data = sample_ij, factorize = FALSE)
  afex_params_ij <- modelbased::estimate_contrasts(
    afex_mod_ij, contrast = "gw", fixed = "cov", adjust = "none") |> as.data.frame()
  names(afex_params_ij) <- paste0("afex_params_", names(afex_params_ij))
  
  
  x1_y1 = sample_ij[sample_ij$gw == 0, ]
  x2_y2 = sample_ij[sample_ij$gw == 1, ]
  
  
  
  rob_trim_output <- 
    tryCatch(
      expr = {
        WRS::Dancova(x1_y1$cov, x1_y1$y, x2_y2$cov, x2_y2$y, tr = 0.2, plotit = FALSE)$output
      }, 
      error = function(e) paste0(e)
    )
    
  if(is.character(rob_trim_output)){
    
    rob_trim_mod_ij <- data.frame(X = NA, n = NA, DIF = NA, TEST = NA, se = NA, ci.low = NA, ci.hi = NA, p.value = NA, p.adjust = NA, 
                                  error = rob_trim_output)
    
  } else {
  
  rob_trim_output_index <- floor(median(1:nrow(rob_trim_output)))
  
  rob_trim_mod_ij <- rob_trim_output[rob_trim_output_index, ] |> t() |> 
    as.data.frame() |> 
    dplyr::mutate(error = NA_character_)
  
  }
  names(rob_trim_mod_ij) <- paste0("rob_trim_params_", names(rob_trim_mod_ij))
  
  
  
  
  
  
  rob_boot_output <- 
    tryCatch(
      expr = {
        WRS::Dancovapb(x1_y1$cov, x1_y1$y, x2_y2$cov, x2_y2$y, tr = 0, plotit = FALSE, SEED = FALSE)$output
      }, 
      error = function(e) paste0(e)
    )
  
  if(is.character(rob_boot_output)){
    
    rob_boot_mod_ij <- data.frame(X = NA, n = NA, DIF = NA, TEST = NA, se = NA, ci.low = NA, ci.hi = NA, p.value = NA, p.adjust = NA, 
                                  error = rob_boot_output)
    
  } else {
    
    rob_boot_output_index <- floor(median(1:nrow(rob_boot_output)))
    
    rob_boot_mod_ij <- rob_boot_output[rob_boot_output_index, ] |> t() |> 
      as.data.frame() |> 
      dplyr::mutate(error = NA_character_)
    
  }
  names(rob_boot_mod_ij) <- paste0("rob_boot_params_", names(rob_boot_mod_ij))
  

  
  
  
  rob_trimboot_output <- 
    tryCatch(
      expr = {
        WRS::Dancovapb(x1_y1$cov, x1_y1$y, x2_y2$cov, x2_y2$y, tr = 0.2, plotit = FALSE, SEED = FALSE)$output
      }, 
      error = function(e) paste0(e)
    )
  
  if(is.character(rob_trimboot_output)){
    
    rob_trimboot_mod_ij <- data.frame(X = NA, n = NA, DIF = NA, TEST = NA, se = NA, ci.low = NA, ci.hi = NA, p.value = NA, p.adjust = NA, 
                                  error = rob_trimboot_output)
    
  } else {
    
    rob_trimboot_output_index <- floor(median(1:nrow(rob_trimboot_output)))
    
    rob_trimboot_mod_ij <- rob_trimboot_output[rob_trimboot_output_index, ] |> t() |> 
      as.data.frame() |> 
      dplyr::mutate(error = NA_character_)
    
  }
  names(rob_trimboot_mod_ij) <- paste0("rob_trimboot_params_", names(rob_trimboot_mod_ij))
  
  
  
  
  dplyr::bind_cols(
    afex_summary_ij, afex_params_ij, 
    rob_trim_mod_ij, rob_boot_mod_ij, rob_trimboot_mod_ij
  )
}


## Iterate -----------------------------------------------------------------

start_time = Sys.time()

set.seed(global_seed)

cat_2w1num_unique_populations <- cat_2w1num_combinations |> 
  dplyr::filter(!duplicated(population_id))

cat_2w1num_n_pop <- nrow(cat_2w1num_unique_populations)


for(i in 1:cat_2w1num_n_pop){
  
  show_progress(i = i, n_pop = cat_2w1num_n_pop)
  
  pop_row_i <- cat_2w1num_unique_populations[i, ]
  pop_i <- generate_cat2w1num(
    n = n_population, b = pop_row_i$b,  
    min_bound = pop_row_i$min_bound, max_bound = pop_row_i$max_bound, 
    vp_fun = text_as_fun(pop_row_i$vp_fun), 
    mod_fun = text_as_fun(pop_row_i$mod_fun),
    vp_val = pop_row_i$vp_val, 
    sigma_shift = pop_row_i$sigma_shift, 
    nu = pop_row_i$nu, tau = pop_row_i$tau
  )
  
  pop_i_wide <- pop_i |> tidyr::pivot_wider(id_cols = c(id, cov), names_from = gw, values_from = y)
  
  grid_i <- cat_2w1num_combinations |> dplyr::filter(population_id == pop_row_i$population_id)
  
  ## start furrr level: sample fit and save
  
  tictoc::tic()
  future::plan(multisession)
  
  furrr_seed <- ifelse((start_iter == 1 & n_iter == 500), global_seed-1, global_seed)
  
  sim_results = furrr::future_pmap(
    .l = grid_i |> dplyr::select(n_rm),
    .f = cat_2w1num_sfs,
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

