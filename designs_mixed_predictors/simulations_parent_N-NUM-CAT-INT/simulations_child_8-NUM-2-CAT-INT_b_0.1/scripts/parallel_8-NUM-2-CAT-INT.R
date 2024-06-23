source("scripts/parallel_setup_script.R")

start_pop_i = 32

# N-NUM-CAT-INT -----------------------------------------------------------

n_1num_1cat_int <- read.csv("data/simulation_grids/N-1NUM-1CAT-INT.csv")
n_1num_2cat_int <- read.csv("data/simulation_grids/N-1NUM-2CAT-INT.csv")

n_3num_1cat_int <- read.csv("data/simulation_grids/N-3NUM-1CAT-INT.csv")
n_3num_2cat_int <- read.csv("data/simulation_grids/N-3NUM-2CAT-INT.csv")

n_8num_1cat_int <- read.csv("data/simulation_grids/N-8NUM-1CAT-INT.csv")
n_8num_2cat_int <- read.csv("data/simulation_grids/N-8NUM-2CAT-INT.csv")

n_num_cat_int <- rbind(
  n_1num_1cat_int |> dplyr::mutate(n_int = 1), 
  n_1num_2cat_int |> dplyr::mutate(n_int = 2), 
  n_3num_1cat_int |> dplyr::mutate(n_int = 1), 
  n_3num_2cat_int |> dplyr::mutate(n_int = 2), 
  n_8num_1cat_int |> dplyr::mutate(n_int = 1), 
  n_8num_2cat_int |> dplyr::mutate(n_int = 2)
)

n_num_cat_int_combinations <- dplyr::full_join(
  tidyr::expand_grid(
    design = "N-NUM-CAT-INT", 
    combo_id = n_num_cat_int$combo_id,
    b, n_bw, n_ratio = n_ratio_bw, 
    #n_int = c(1,2),
    n_cont_preds = c(1,3,8), 
    n_cat_preds = c(1,2),
    iter = start_iter:n_iter
  ), 
  n_num_cat_int |> dplyr::select(-design),
  by = c("combo_id", "n_cat_preds", "n_cont_preds")
) |> 
  dplyr::filter(
    n_cont_preds == 8, 
    n_cat_preds == 2, 
    b == 0.1
  ) |> 
  dplyr::mutate(
    population_id = paste0(design, "_", combo_id, "_b_", b, 
                           "_n_cont_preds_", n_cont_preds, "_n_cat_preds_", n_cat_preds, 
                           "_n_int_", n_int)
  )


generate_n_num_cat_int <- function(
    n = 100000,
    b,
    n_cont_preds,
    n_cat_preds,
    n_int,
    min_bound, 
    max_bound, 
    nu, 
    tau, 
    vp_val, 
    sigma_shift, 
    n_het_preds, 
    mod_fun,
    vp_fun
    
){
  
  sigma = 1
  r = 0.212
  r_cat_shift = 0.042
  
  cat_pred_names = paste0("x", 1:n_cat_preds, "_cat")
  
  # generate correlated predictors 
  set.seed(global_seed) 
  norm_df <- faux::rnorm_multi(n = n, r = r, varnames = paste0("x", 1:n_cont_preds))  
  
  set.seed(global_seed)
  norm_df[[cat_pred_names[1]]] <- faux::rnorm_pre(norm_df,  r = r + r_cat_shift)
  
  if(n_cat_preds == 2){
    norm_df[[cat_pred_names[2]]] <- faux::rnorm_pre(norm_df, r = r + r_cat_shift)
    norm_df <- norm_df |> dplyr::arrange(by = !!sym(cat_pred_names[2]))
    norm_df[[cat_pred_names[2]]] <- rep(c(0, 1), each = n/2)
  }
  
  norm_df <- norm_df |> dplyr::arrange(by = !!sym(cat_pred_names[1]))
  norm_df[[cat_pred_names[1]]] <- rep(c(0, 1), each = n/2)
  
  # generate errors 
  set.seed(21023) 
  e <- gamlss.dist::rSHASHo(
    n = n,
    mu = 0,
    nu = nu,
    tau = tau
  )
  
  # modify error 
  if(n_het_preds == 1){
    
    x1_mod <- mod_fun(norm_df$x1)
    x_mult <- (x1_mod)^(1/1)
    m_cont <- vp_fun(x = x_mult, vp_val)
    
    x1_cat_sd = ifelse(norm_df$x1_cat == 1, sigma + sigma_shift, sigma)
    m_cat <- x1_cat_sd
    
    m_total <- m_cat*m_cont
    
  } else if(n_het_preds == 2){
    
    x1_mod <- mod_fun(norm_df$x1)
    x_mult <- (x1_mod)^(1/1)
    m_cont <- vp_fun(x = x_mult, vp_val)
    
    x1_cat_sd = ifelse(norm_df$x1_cat == 1, sigma + sigma_shift, sigma)
    m_cat <- x1_cat_sd
    
    m_total <- m_cat*m_cont
    
  } else if(n_het_preds == 0){
    m_total = 1
  }
  
  #m_total = 1
  e_mod = e*m_total
  e_mod <- squish(e_mod, min_bound, max_bound)
  
  # generate y 
  if(n_int == 1){
    
    norm_df$y <- 
      rowSums(norm_df[, 1:(n_cont_preds+n_cat_preds)]*b) + 
      b*(norm_df$x1*norm_df$x1_cat) + e_mod
    
  } else if(n_int == 2){
    
    norm_df$y <- 
      rowSums(norm_df[, 1:(n_cont_preds+n_cat_preds)]*b) + 
      b*(norm_df$x1*norm_df$x1_cat) + 
      b*(norm_df$x1*norm_df$x2_cat) + e_mod
    
  }
  
  return(norm_df)
}


## sample fit and save --------------------------------------------------------

n_num_cat_int_sfs <- function(
   # pop_i_x1_cat_0, 
   # pop_i_x1_cat_1, 
    n_cat_preds, n_int, n_bw, n_ratio){
  
  ## create a sample   
  ns <- sample_n_from_ratio_bw2(n_bw = n_bw, n_ratio = n_ratio)
  
  if(n_cat_preds == 1){
    
    sample_ij <- dplyr::bind_rows(
      pop_i_x1_cat_0 |> dplyr::slice_sample(n = ns[["n1"]]), 
      pop_i_x1_cat_1 |> dplyr::slice_sample(n = ns[["n2"]]) 
    )
    
  } else if(n_cat_preds == 2){
    
    weight_col <- paste0("n_ratio_w_", n_ratio)
    
    sample_ij <- dplyr::bind_rows(
      pop_i_x1_cat_0 |> dplyr::slice_sample(n = ns[["n1"]], weight_by = (!!sym(weight_col))), 
      pop_i_x1_cat_1 |> dplyr::slice_sample(n = ns[["n2"]], weight_by = (!!sym(weight_col))) 
    )
    
    n_cell_count <- sample_ij |> dplyr::count(x1_cat, x2_cat) |> dplyr::pull(n)
    
    # must sample at least 2 for each cell and all groups must be represented
    # in the cell grid 
    while (any(n_cell_count < 2) | length(n_cell_count) < 4) {
      
      sample_ij <- dplyr::bind_rows(
        pop_i_x1_cat_0 |> dplyr::slice_sample(n = ns[["n1"]], weight_by = (!!sym(weight_col))), 
        pop_i_x1_cat_1 |> dplyr::slice_sample(n = ns[["n2"]], weight_by = (!!sym(weight_col))) 
      )
      
      n_cell_count <- sample_ij |> dplyr::count(x1_cat, x2_cat) |> dplyr::pull(n)
      
    }
    
  }
  
  sample_ij <- sample_ij |> 
    dplyr::select(-dplyr::contains("ratio"))
  
  
  if(n_int == 1){
    
    lm_formula = paste0(
      "y ~ ", 
      paste0(names(sample_ij[, -ncol(sample_ij)]), collapse = " + "), 
      " + x1:x1_cat"
    )
    
  } else if(n_int == 2){
    
    lm_formula = paste0(
      "y ~ ", 
      paste0(names(sample_ij[, -ncol(sample_ij)]), collapse = " + "), 
      " + x1:x1_cat + x1:x2_cat"
    )
    
  }
  
  # Fit models and export: 
  fit_and_summarise_mods_bw(mod_formula = lm_formula, sample_ij = sample_ij)
  
}


## Iterate -----------------------------------------------------------------

  
set.seed(global_seed)

n_num_cat_int_unique_populations <- n_num_cat_int_combinations |> 
  dplyr::filter(!duplicated(population_id))

n_num_cat_int_n_pop <- nrow(n_num_cat_int_unique_populations)


for(i in start_pop_i:n_num_cat_int_n_pop){
  
  pop_row_i <- n_num_cat_int_unique_populations[i, ]
  
  grid_i <- n_num_cat_int_combinations |> dplyr::filter(population_id == pop_row_i$population_id)
  
  n_ratios <- sort(unique(grid_i$n_ratio))
  
  pop_i <- generate_n_num_cat_int(
    n_int = pop_row_i$n_int,
    sigma_shift = pop_row_i$sigma_shift,
    n_cat_preds = pop_row_i$n_cat_preds, 
    n_cont_preds = pop_row_i$n_cont_preds,
    n = n_population, b = pop_row_i$b, 
    min_bound = pop_row_i$min_bound, max_bound = pop_row_i$max_bound, 
    vp_val = pop_row_i$vp_val, n_het_preds = pop_row_i$n_het_pred, 
    mod_fun = text_as_fun(pop_row_i$mod_fun), 
    vp_fun = text_as_fun(pop_row_i$vp_fun), 
    nu = pop_row_i$nu, tau = pop_row_i$tau
  )
  
  if(pop_row_i$n_cat_preds == 2){
    
    pop_i[[paste0("n_ratio_w_", n_ratios[1])]] <- ifelse(pop_i$x2_cat == "1" , n_ratios[1], 1)
    pop_i[[paste0("n_ratio_w_", n_ratios[2])]] <- ifelse(pop_i$x2_cat == "1" , n_ratios[2], 1)
    pop_i[[paste0("n_ratio_w_", n_ratios[3])]] <- ifelse(pop_i$x2_cat == "1" , n_ratios[3], 1)
    
  }
  
  pop_i_x1_cat_0 <- pop_i |> dplyr::filter(x1_cat == 0)
  pop_i_x1_cat_1 <- pop_i |> dplyr::filter(x1_cat == 1)
  
  grid_i <- n_num_cat_int_combinations |> dplyr::filter(population_id == pop_row_i$population_id)
  
  ## start furrr level: sample fit and save
  
  tictoc::tic()
  future::plan(multisession)
  
  furrr_seed <- ifelse((start_iter == 1 & n_iter == 500), global_seed-1, global_seed)
  
  sim_results = furrr::future_pmap(
    .l = grid_i |> dplyr::select(n_cat_preds, n_int, n_bw, n_ratio), 
    .f = n_num_cat_int_sfs,
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
  rm(pop_i_x1_cat_0)
  rm(pop_i_x1_cat_1)
  
}

rm(pop_row_i)
rm(pop_i)
rm(grid_i)
gc()


