source("scripts/parallel_setup_script.R")

start_pop_i = 1

# CAT-2W-2B ---------------------------------------------------------------

cat_2w2b <- readr::read_csv("data/simulation_grids/CAT-2W-2B.csv")

cat_2w2b_combinations <- dplyr::full_join(
  tidyr::expand_grid(
    design = "CAT-2W-2B", 
    combo_id = cat_2w2b$combo_id,
    b, n_rm, 
    n_ratio = n_ratio_rm, 
    iter = start_iter:n_iter
  ), 
  cat_2w2b |> dplyr::select(-design),
  by = c("combo_id")
) |> 
  dplyr::mutate(
    population_id = paste0(design, "_", combo_id, "_b_", b)
  )

generate_cat2w2b <- function(
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
    within = list(w_group1 = c("gw1", "gw2")),
    between = list(between_group = c("gb1", "gb2")),
    plot = FALSE
  )
  
  transformed_df <- norm_df |> 
    dplyr::mutate(
      gw1 = dplyr::case_when(
        between_group == "gb1" ~ transform_norm(raw_var = gw1, mu = 0, nu = nu, tau = tau,
                                                sigma = sigma),
        between_group == "gb2" ~ transform_norm(raw_var = gw1, mu = 0, nu = nu, tau = tau,
                                                sigma = sigma + sigma_shift_1) + b + b_shift_1
      ),
      gw2 = dplyr::case_when(
        between_group == "gb1" ~ transform_norm(raw_var = gw2, mu = 0, nu = nu, tau = tau,
                                                sigma = sigma) + b,
        between_group == "gb2" ~ transform_norm(raw_var = gw2, mu = 0, nu = nu, tau = tau,
                                                sigma = sigma + sigma_shift_2) + 3*b + b_shift_2
      ),
    )
  
  return(transformed_df)
}

tsplit_wrapper <- function(J, K, blist, tr){
  warn_increase <- 0
  warn_message <- NULL
  
  tsplit_result <- withCallingHandlers(
    
    WRS::tsplit(J, K, blist, tr),
    
    warning = function(w) {
      warn_increase <<- 1
      warn_message <<- paste0(w)
    }
  )
  
  tsplit_result$warn_increase <- warn_increase
  tsplit_result$warn_message <- warn_message
  
  return(tsplit_result)
}

tsplit_error_wrapper <- function(J, K, blist, tr){
  tryCatch(
    tsplit_wrapper(J, K, blist, tr) |> suppressWarnings(),
    error = function(e) {
      paste0(e)
    }
  )
}

bwtrimbt_mod <- function (J, K, x, tr = 0.2, JK = J * K, grp = c(1:JK), nboot = 599, SEED = FALSE) 
{
  if (SEED)
    set.seed(2)
  if (is.data.frame(x) || is.matrix(x)) {
    y <- list()
    ik = 0
    il = c(1:K) - K
    for (j in 1:J) {
      il = il + K
      zz = x[, il]
      zz = elimna(zz)
      for (k in 1:K) {
        ik = ik + 1
        y[[ik]] = zz[, k]
      }
    }
    x <- y
  }
  JK <- J * K
  data <- list()
  xcen <- list()
  for (j in 1:length(x)) {
    data[[j]] <- x[[grp[j]]]
    xcen[[j]] <- data[[j]] - mean(data[[j]], tr)
  }
  x <- data
  # set.seed(2) # Why would you do this to me.  
  nvec <- NA
  jp <- 1 - K
  for (j in 1:J) {
    jp <- jp + K
    nvec[j] <- length(x[[j]])
  }
  
  blist <- list()
  
  warn_count = 0
  warn_vec <- NULL
  error_count = 0
  error_vec <- NULL
  #print("Taking bootstrap samples. Please wait.") #No thanks. 
  testmat <- matrix(NA, ncol = 3, nrow = nboot)
  for (iboot in 1:nboot) {
    iv <- 0
    for (j in 1:J) {
      temp <- sample(nvec[j], replace = T)
      for (k in 1:K) {
        iv <- iv + 1
        tempx <- xcen[[iv]]
        blist[[iv]] <- tempx[temp]
      }
    }
    btest <- tsplit_error_wrapper(J, K, blist, tr)
    
    if(is.character(btest)){
      
      testmat[iboot, 1] <- NA
      testmat[iboot, 2] <- NA
      testmat[iboot, 3] <- NA
      
      error_count = error_count + 1
      error_vec <- c(error_vec, btest)
      
    } else {
      
      warn_count = c(warn_count, btest$warn_increase)
      warn_vec = c(warn_vec, btest$warn_message)
      
      testmat[iboot, 1] <- btest$Qa
      testmat[iboot, 2] <- btest$Qb
      testmat[iboot, 3] <- btest$Qab
      
    }
  }
  
  test = WRS::tsplit(J, K, x, tr = tr)
  pbA = mean(test$Qa[1] < testmat[, 1])
  pbB = mean(test$Qb[1] < testmat[, 2])
  pbAB = mean(test$Qab[1] < testmat[, 3])
  
  list(
    p.value.A = pbA, p.value.B = pbB, p.value.AB = pbAB, 
    warn_count = sum(warn_count), 
    warn_message = ifelse(is.null(warn_vec), NA_character_, paste0(unique(warn_vec), sep = ";")), 
    error_count = error_count, 
    error_message = ifelse(is.null(error_vec), NA_character_, paste0(unique(error_vec), sep = ";"))
  )
}

## sample fit and save -------------------------------------------------------

cat_2w2b_sfs <- function(
    #pop_i_gb1, 
    #pop_i_gb2 , 
    n_rm, 
    n_ratio
    ){
  
  ## calculate sample split 
  n1n2 <- sample_n_from_ratio_rm(n_rm, n_ratio)
  # n1n2 <- sample_n_from_ratio_rm(18, 2)
  
  ## create a sample   
  sample_ij <- dplyr::bind_rows(
    pop_i_gb1 |> dplyr::slice_sample(n = n1n2[["n1"]]), 
    pop_i_gb2 |> dplyr::slice_sample(n = n1n2[["n2"]]) 
  )
  
  
  ## fit ols model 
  
  sample_ij_long <- sample_ij |> 
    tidyr::pivot_longer(-c(id, between_group), names_to = "gw", values_to = "y") 
  
  afex_summary_ij <- afex_nc_gg_summary("y ~ gw*between_group + (gw|id)", data = sample_ij_long)
  
  ## robust mods 
  
  wrs_list <- list(
    `1a` = dplyr::filter(sample_ij_long, between_group == "gb1", gw == "gw1")$y,
    `1b` = dplyr::filter(sample_ij_long, between_group == "gb1", gw == "gw2")$y,
    `2a` = dplyr::filter(sample_ij_long, between_group == "gb2", gw == "gw1")$y,
    `2b` = dplyr::filter(sample_ij_long, between_group == "gb2", gw == "gw2")$y
  )
  
  if(any(n1n2 < 12)){nboot = 1000
  } else{nboot = 599}
  
  rob_trim_mod_ij <- WRS::bwtrim(J = 2, K = 2, data = wrs_list, tr = 0.2)
  rob_trim_summary_ij <- tibble::tibble(
    statistic = unlist(rob_trim_mod_ij[c(1,3,5)]),
    p = unlist(rob_trim_mod_ij[c(2,4,6)]))
  names(rob_trim_summary_ij) <- paste0("rob_trim_", names(rob_trim_summary_ij))
  
  rob_boot_mod_ij <-  bwtrimbt_mod(J = 2, K = 2,  x = wrs_list, tr = 0, nboot = nboot) |> suppressWarnings()
  rob_boot_summary_ij <- tibble::tibble(
    rob_boot_p = unlist(rob_boot_mod_ij[-c(4,5,6,7)]), 
    rob_boot_warn_count = unlist(rob_boot_mod_ij[4]), 
    rob_boot_warns = unlist(rob_boot_mod_ij[5]), 
    rob_boot_error_count = unlist(rob_boot_mod_ij[6]),
    rob_boot_errors = unlist(rob_boot_mod_ij[7])
  ) 
  
  rob_trimboot_mod_ij <-  bwtrimbt_mod(J = 2, K = 2,  x = wrs_list, tr = 0.2, nboot = nboot) |> suppressWarnings()
  rob_trimboot_summary_ij <- tibble::tibble(
    rob_trimboot_p = unlist(rob_trimboot_mod_ij[-c(4,5,6,7)]), 
    rob_trimboot_warn_count = unlist(rob_trimboot_mod_ij[4]), 
    rob_trimboot_warns = unlist(rob_trimboot_mod_ij[5]), 
    rob_trimboot_error_count = unlist(rob_trimboot_mod_ij[6]),
    rob_trimboot_errors = unlist(rob_trimboot_mod_ij[7])
  )
  
  afex_int_b <- WRS::bwmcp(2,2, wrs_list, tr = 0, nboot = 2)$Fac.AB[2]
  rob_trim_int_b <- WRS::bwmcp(2,2, wrs_list, tr = 0.2, nboot = 2)$Fac.AB[2]
  
  dplyr::bind_cols(
    afex_summary_ij, 
    rob_trim_summary_ij, 
    rob_boot_summary_ij, 
    rob_trimboot_summary_ij, 
    afex_int_b = afex_int_b, rob_trim_int_b = rob_trim_int_b
  )
  
}



## Iterate ----------------------------------------------------------------

set.seed(global_seed)

cat_2w2b_unique_populations <- cat_2w2b_combinations |> 
  dplyr::filter(!duplicated(population_id))

cat_2w2b_n_pop <- nrow(cat_2w2b_unique_populations)

for(i in start_pop_i:cat_2w2b_n_pop){
  
  show_progress(i = i, n_pop = cat_2w2b_n_pop)
  
  pop_row_i <- cat_2w2b_unique_populations[i, ]
  pop_i <- generate_cat2w2b(
    n = n_population, b = pop_row_i$b, 
    sigma_shift_1 = pop_row_i$sigma_shift_1, sigma_shift_2 = pop_row_i$sigma_shift_2, 
    nu = pop_row_i$nu, tau = pop_row_i$tau
  )
  
  pop_i_gb1 <- pop_i |> dplyr::filter(between_group == "gb1") 
  pop_i_gb2 <- pop_i |> dplyr::filter(between_group == "gb2") 
  
  grid_i <- cat_2w2b_combinations |> dplyr::filter(population_id == pop_row_i$population_id)
  
  ## start furrr level: sample fit and save
  
  tictoc::tic()
  future::plan(multisession)
  
  furrr_seed <- ifelse((start_iter == 1 & n_iter == 500), global_seed-1, global_seed)
  
  sim_results = furrr::future_pmap(
    .l = grid_i |> dplyr::select(n_rm, n_ratio), 
    .f = cat_2w2b_sfs,
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
  
  rm(pop_i_gb1)
  rm(pop_i_gb2)
  
}

rm(pop_row_i)
rm(pop_i)
rm(grid_i)
gc()



