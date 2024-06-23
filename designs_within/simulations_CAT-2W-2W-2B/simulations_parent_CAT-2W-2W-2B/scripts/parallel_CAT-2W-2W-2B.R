source("scripts/parallel_setup_script.R")

# CAT-2W-2W-2B ------------------------------------------------------------

cat_2w2w2b <- readr::read_csv("data/simulation_grids/CAT-2W-2W-2B.csv")

cat_2w2w2b_combinations <- dplyr::full_join(
  tidyr::expand_grid(
    design = "CAT-2W-2W-2B", 
    combo_id = cat_2w2w2b$combo_id,
    b, n_rm, 
    n_ratio = n_ratio_rm, 
    iter = start_iter:n_iter
  ), 
  cat_2w2w2b |> dplyr::select(-design),
  by = c("combo_id")
) |> 
  dplyr::mutate(
    population_id = paste0(design, "_", combo_id, "_b_", b)
  )


generate_cat2w2w2b <- function(
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
    within = list(
      factor_1 = c("1", "2"), factor_2 = c("a", "b")
    ), 
    between = list(gb1 = c(".", "..")),
    plot = FALSE
  )
  
  transformed_df <- norm_df |> 
    dplyr::mutate(
      `1_a` = dplyr::case_when(
        gb1 ==  "." ~ transform_norm(raw_var = `1_a`, mu = 0, nu = nu, tau = tau, sigma = sigma),
        gb1 == ".." ~ transform_norm(raw_var = `1_a`, mu = 0, nu = nu, tau = tau, sigma = sigma) + b,
      ), 
      `1_b` = dplyr::case_when(
        gb1 ==  "." ~ transform_norm(raw_var = `1_b`, mu = 0, nu = nu, tau = tau, sigma = sigma) + b,
        gb1 == ".." ~ transform_norm(raw_var = `1_b`, mu = 0, nu = nu, tau = tau, sigma = sigma + sigma_shift_1) + 3*b + b_shift_1,
      ), 
      `2_a` = dplyr::case_when(
        gb1 ==  "." ~ transform_norm(raw_var = `2_a`, mu = 0, nu = nu, tau = tau, sigma = sigma) + b,
        gb1 == ".." ~ transform_norm(raw_var = `2_a`, mu = 0, nu = nu, tau = tau, sigma = sigma) + 3*b,
      ), 
      `2_b` = dplyr::case_when(
        gb1 ==  "." ~ transform_norm(raw_var = `2_b`, mu = 0, nu = nu, tau = tau, sigma = sigma) + 3*b,
        gb1 == ".." ~ transform_norm(raw_var = `2_b`, mu = 0, nu = nu, tau = tau, sigma = sigma + sigma_shift_2) + 7*b + b_shift_2,
      )
    )
  
  return(transformed_df)
}

bwwtrim_wrapper <- function(J, K, L, bsam, tr){
  warn_increase <- 0
  warn_message <- NULL
  
  bwwtrim_result <- withCallingHandlers(
    
    WRS::bwwtrim(J, K, L, bsam, tr = tr),
    
    warning = function(w) {
      warn_increase <<- 1
      warn_message <<- paste0(w)
    }
  )
  
  bwwtrim_result$warn_increase <- warn_increase
  bwwtrim_result$warn_message <- warn_message
  
  return(bwwtrim_result)
}

bwwtrim_error_wrapper <- function(J, K, L, bsam, tr){
  tryCatch(
    bwwtrim_wrapper(J, K, L, bsam, tr) |> suppressWarnings(),
    error = function(e) {
      paste0(e)
    }
  )
}

bwwtrimbt_mod <- function (J, K, L, x, tr = 0.2, JKL = J * K * L, con = 0, alpha = 0.05, 
                           grp = c(1:JKL), nboot = 599, SEED = FALSE, ...) 
{
  if (is.data.frame(x)) 
    data = as.matrix(x)
  if (is.matrix(x)) {
    y <- list()
    for (j in 1:ncol(x)) y[[j]] <- x[, j]
    x <- y
  }
  conM = con3way(J, K, L)
  p <- J * K * L
  if (p > length(x)) 
    stop("JKL is less than the Number of groups")
  JK = J * K
  KL = K * L
  v <- matrix(0, p, p)
  data <- list()
  xx = list()
  for (j in 1:length(x)) {
    data[[j]] <- x[[grp[j]]]
    xx[[j]] = x[[grp[j]]]
    data[[j]] = data[[j]] - mean(data[[j]], tr = tr)
  }
  x <- data
  test = bwwtrim(J, K, L, xx, tr = tr)
  if (SEED) 
    set.seed(2)
  bsam = list()
  bdat = list()
  
  warn_count = 0
  warn_vec <- NULL
  error_count = 0
  error_vec <- NULL
  
  aboot = NA
  bboot = NA
  cboot = NA
  abboot = NA
  acboot = NA
  bcboot = NA
  abcboot = NA
  for (ib in 1:nboot) {
    ilow <- 1 - KL
    iup = 0
    for (j in 1:J) {
      ilow <- ilow + KL
      iup = iup + KL
      nv = length(x[[ilow]])
      bdat[[j]] = sample(nv, size = nv, replace = TRUE)
      for (k in ilow:iup) {
        bsam[[k]] = x[[k]][bdat[[j]]]
      }
    }
    temp = bwwtrim_error_wrapper(J, K, L, bsam, tr = tr)
    
    if(is.character(temp)){
      
      aboot[ib] = NA
      bboot[ib] = NA
      cboot[ib] = NA
      acboot[ib] = NA
      bcboot[ib] = NA
      abboot[ib] = NA
      abcboot[ib] = NA
      
      error_count = 0
      error_vec <- c(error_vec, temp)
      
    } else {
      
      warn_count = c(warn_count, temp$warn_increase)
      warn_vec = c(warn_vec, temp$warn_message)
      
      aboot[ib] = temp$Qa
      bboot[ib] = temp$Qb
      cboot[ib] = temp$Qc
      acboot[ib] = temp$Qac
      bcboot[ib] = temp$Qbc
      abboot[ib] = temp$Qab
      abcboot[ib] = temp$Qabc
      
    }
    
  }
  
  pbA = NA
  pbB = NA
  pbC = NA
  pbAB = NA
  pbAC = NA
  pbBC = NA
  pbABC = NA
  pbA = mean(test$Qa[1, 1] < aboot, na.rm = TRUE)
  pbB = mean(test$Qb[1, 1] < bboot, na.rm = TRUE)
  pbC = mean(test$Qc[1, 1] < cboot, na.rm = TRUE)
  pbAB = mean(test$Qab[1, 1] < abboot, na.rm = TRUE)
  pbAC = mean(test$Qac[1, 1] < acboot, na.rm = TRUE)
  pbBC = mean(test$Qbc[1, 1] < bcboot, na.rm = TRUE)
  pbABC = mean(test$Qabc[1, 1] < abcboot, na.rm = TRUE)
  
  list(p.value.A = pbA, p.value.B = pbB, p.value.C = pbC, 
       p.value.AB = pbAB, p.value.AC = pbAC, p.value.BC = pbBC, 
       p.value.ABC = pbABC, 
       warn_count = sum(warn_count), 
       warn_message = ifelse(is.null(warn_vec), NA_character_, paste0(unique(warn_vec), sep = ";")), 
       error_count = error_count, 
       error_message = ifelse(is.null(error_vec), NA_character_, paste0(unique(error_vec), sep = ";"))
  )
}


## sample fit and save
cat_2w2w2b_sfs <- function(
   # pop_i_., pop_i_.., 
    n_rm, n_ratio){
  
  ## calculate sample split 
  n1n2 <- sample_n_from_ratio_rm(n_rm, n_ratio)
  #n1n2 <- sample_n_from_ratio_rm(18, 2)
  
  ## create a sample   
  sample_ij <- dplyr::bind_rows(
    pop_i_. |> dplyr::slice_sample(n = n1n2[["n1"]]), 
    pop_i_.. |> dplyr::slice_sample(n = n1n2[["n2"]]) 
  )
  
  ## fit ols model 
  
  sample_ij_long <- sample_ij |> 
    tidyr::pivot_longer(-c(id, gb1), names_to = "gw", values_to = "y") |> 
    dplyr::mutate(
      gw1 = ifelse(stringr::str_detect(gw, "1"), "1", "2"), 
      gw2 = ifelse(stringr::str_detect(gw, "a"), "a", "b")
    ) 
  
  afex_summary_ij <-  afex_nc_gg_summary("y ~ gb1*gw1*gw2 + (gw1*gw2|id)", data = sample_ij_long) |> 
    dplyr::arrange(afex_nc_term)
  
  wrs_df <- list(
    dplyr::filter(sample_ij_long, gb1 == ".", gw1 == "1", gw2 == "a")$y, 
    dplyr::filter(sample_ij_long, gb1 == ".", gw1 == "1", gw2 == "b")$y, 
    dplyr::filter(sample_ij_long, gb1 == ".", gw1 == "2", gw2 == "a")$y, 
    dplyr::filter(sample_ij_long, gb1 == ".", gw1 == "2", gw2 == "b")$y,
    
    dplyr::filter(sample_ij_long, gb1 == "..", gw1 == "1", gw2 == "a")$y, 
    dplyr::filter(sample_ij_long, gb1 == "..", gw1 == "1", gw2 == "b")$y, 
    dplyr::filter(sample_ij_long, gb1 == "..", gw1 == "2", gw2 == "a")$y, 
    dplyr::filter(sample_ij_long, gb1 == "..", gw1 == "2", gw2 == "b")$y
  )
  
  if(any(n1n2 < 12)){nboot = 1000
  } else{nboot = 599}
  
  rob_trim_mod_ij <- WRS::bwwtrim(2,2,2, wrs_df, tr = 0.2)
  rob_trim_summary_ij <- tibble::tibble(
    rob_trim_term = names(rob_trim_mod_ij[grep(pattern = "p.value",  names(rob_trim_mod_ij), invert = TRUE)]),
    rob_trim_statistic = unlist(rob_trim_mod_ij[grep(pattern = "p.value",  names(rob_trim_mod_ij), invert = TRUE)]), 
    rob_trim_p = unlist(rob_trim_mod_ij[grep(pattern = "p.value",  names(rob_trim_mod_ij), invert = FALSE)])
  ) |> dplyr::arrange(rob_trim_term)
  
  rob_boot_mod_ij <- bwwtrimbt_mod(2,2,2, wrs_df, tr = 0, SEED = FALSE, nboot = nboot) |> suppressWarnings()
  rob_boot_summary_ij <- tibble::tibble(
    rob_boot_term = names(rob_boot_mod_ij[grep(pattern = "p.value",  names(rob_boot_mod_ij), invert = FALSE)]),
    rob_boot_p = unlist(rob_boot_mod_ij[grep(pattern = "p.value",  names(rob_boot_mod_ij), invert = FALSE)]), 
    rob_boot_warn_count = unlist(rob_boot_mod_ij[grep(pattern = "warn_count",  names(rob_boot_mod_ij), invert = FALSE)]), 
    rob_boot_warn_message = unlist(rob_boot_mod_ij[grep(pattern = "warn_message",  names(rob_boot_mod_ij), invert = FALSE)]),
    rob_boot_error_count = unlist(rob_boot_mod_ij[grep(pattern = "error_count",  names(rob_boot_mod_ij), invert = FALSE)]), 
    rob_boot_error_message = unlist(rob_boot_mod_ij[grep(pattern = "error_message",  names(rob_boot_mod_ij), invert = FALSE)])
  ) |> dplyr::arrange(rob_boot_term)
  
  rob_trimboot_mod_ij <- bwwtrimbt_mod(2,2,2, wrs_df, tr = 0.2, SEED = FALSE, nboot = nboot) |> suppressWarnings()
  rob_trimboot_summary_ij <- tibble::tibble(
    rob_trimboot_term = names(rob_trimboot_mod_ij[grep(pattern = "p.value",  names(rob_trimboot_mod_ij), invert = FALSE)]),
    rob_trimboot_p = unlist(rob_trimboot_mod_ij[grep(pattern = "p.value",  names(rob_trimboot_mod_ij), invert = FALSE)]), 
    rob_trimboot_warn_count = unlist(rob_trimboot_mod_ij[grep(pattern = "warn_count",  names(rob_trimboot_mod_ij), invert = FALSE)]), 
    rob_trimboot_warn_message = unlist(rob_trimboot_mod_ij[grep(pattern = "warn_message",  names(rob_trimboot_mod_ij), invert = FALSE)]),
    rob_trimboot_error_count = unlist(rob_trimboot_mod_ij[grep(pattern = "error_count",  names(rob_trimboot_mod_ij), invert = FALSE)]), 
    rob_trimboot_error_message = unlist(rob_trimboot_mod_ij[grep(pattern = "error_message",  names(rob_trimboot_mod_ij), invert = FALSE)])
  ) |> dplyr::arrange(rob_trimboot_term)
  
  afex_int_b <- WRS::bwwmcp(2,2,2, wrs_df, tr = 0, nboot = 2)$Fac.ABC[2]
  rob_trim_int_b <- WRS::bwwmcp(2,2,2, wrs_df, tr = 0.2, nboot = 2)$Fac.ABC[2]
  
  dplyr::bind_cols(
    afex_summary_ij, rob_trim_summary_ij, rob_boot_summary_ij, rob_trimboot_summary_ij, 
    afex_int_b = afex_int_b,rob_trim_int_b = rob_trim_int_b
  )
  
}



## Iterate -----------------------------------------------------------------

set.seed(global_seed)

cat_2w2w2b_unique_populations <- cat_2w2w2b_combinations |> 
  dplyr::filter(!duplicated(population_id))

cat_2w2w2b_n_pop <- nrow(cat_2w2w2b_unique_populations)


for(i in 11:cat_2w2w2b_n_pop){
  
  show_progress(i = i, n_pop = cat_2w2w2b_n_pop)
  
  pop_row_i <- cat_2w2w2b_unique_populations[i, ]
  pop_i <- generate_cat2w2w2b(
    n = n_population, b = pop_row_i$b, 
    sigma_shift_1 = pop_row_i$sigma_shift_1, sigma_shift_2 = pop_row_i$sigma_shift_2, 
    nu = pop_row_i$nu, tau = pop_row_i$tau
  )
  
  pop_i_. <- pop_i |> dplyr::filter(gb1 == ".")
  pop_i_.. <- pop_i |> dplyr::filter(gb1 == "..")
  
  grid_i <- cat_2w2w2b_combinations |> dplyr::filter(population_id == pop_row_i$population_id)
  
  ## start furrr level: sample fit and save
  
  tictoc::tic()
  future::plan(multisession)
  
  furrr_seed <- ifelse((start_iter == 1 & n_iter == 500), global_seed-1, global_seed)
  
  sim_results = furrr::future_pmap(
    .l = grid_i |> dplyr::select(n_rm, n_ratio), 
    .f = cat_2w2w2b_sfs,
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
  
  rm(pop_i_.)
  rm(pop_i_..)
  
}

rm(pop_row_i)
rm(pop_i)
rm(grid_i)
gc()


