source(here::here("scripts/helpers.R"))


# set up files  -----------------------------------------------------------

between_cat <- all_files[grep("designs_between", all_files)]

between_cat_list <- purrr::map(
  .x = between_cat, 
  .f = readRDS
)

names(between_cat_list) <- sub(name_pattern, "", between_cat) |> 
  sub(x = _, replacement = "", pattern = ".rds")


# processing tweaks -------------------------------------------------------

# add nice names to power plots 

between_cat_list <- purrr::map(
  .x = between_cat_list, 
  .f = ~ purrr::map_at(
    .x, 
    .at = "sim_df_power_sum", 
    .f = ~ dplyr::mutate(
      .x, 
      name = factor(name, labels = c("OLS", "Bootstrap", "HC4", "MM/KS", "MM")))
  )
)

# add nice names to power plots 

between_cat_list <- purrr::map(
  .x = between_cat_list, 
  .f = ~ purrr::map_at(
    .x, 
    .at = "sim_df_cover_sum", 
    .f = ~ dplyr::mutate(
      .x, 
      name = factor(name, labels = c("OLS", "Bootstrap", "HC4", "MM/KS", "MM")))
  )
)

# add nice names to error plots 

between_cat_list <- purrr::map(
  .x = between_cat_list, 
  .f = ~ purrr::map_at(
    .x, 
    .at = "sim_df_model_errors_sum", 
    .f = ~ dplyr::mutate(
      .x, 
      name = factor(name, labels = c("OLS", "Bootstrap", "HC4", "MM/KS", "MM")))
  )
)

# generate qmd -----------------------------------------------------------


qmd_lines <- c()

for(i in names(between_cat_list)){
  
  design_name_i <- i
  
  unique_plot_groups_i <- between_cat_list[[i]]$sim_df_overall_sum$plot_group |> unique() |> length()
  
  qmd_lines_i <- paste0(
    
    "## ", design_name_i, "\n\n",
    
    "### Overall model fit", "\n\n",
    
    codechunk(
      settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
      code_string = paste0(
        "power_overall_plot_between_cat(", "between_cat_list[['",i,"']]$sim_df_overall_sum", ")"
      )
    ), "\n\n",

    "### Power/false-positives", "\n\n",

    codechunk(
       settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
      code_string = paste0(
        "power_plot_between_cat(", "between_cat_list[['",i,"']]$sim_df_power_sum", ")"
      )
    ), "\n\n",

    
    "### Bias \n\n",
    
    "#### Mean bias \n\n",

    codechunk(
       settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
      code_string = paste0(
        "bias_mean_plot_between_cat(", "between_cat_list[['",i,"']]$sim_df_bias_sum", ")"
      )
    ), "\n\n",
    
    "#### SD bias \n\n",

    codechunk(
       settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
      code_string = paste0(
        "bias_sd_plot_between_cat(", "between_cat_list[['",i,"']]$sim_df_bias_sum", ")"
      )
    ), "\n\n",
    
    "### CI coverage", "\n\n",

    codechunk(
       settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
      code_string = paste0(
        "cover_plot_between_cat(", "between_cat_list[['",i,"']]$sim_df_cover_sum", ")"
      )
    ), "\n\n",
    
    "### CI bias-shifted coverage", "\n\n",

    codechunk(
       settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
      code_string = paste0(
        "cover_plot_between_cat(", "between_cat_list[['",i,"']]$sim_df_cover_shift_sum", ")"
      )
    ), "\n\n",
    
    "### CI width", "\n\n",

    codechunk(
       settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
      code_string = paste0(
        "ci_width_plot_between_cat(", "between_cat_list[['",i,"']]$sim_df_ci_width_sum", ")"
      )
    ), "\n\n",
    
    
    "### Model errors", "\n\n",

    codechunk(
       settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
      code_string = paste0(
        "model_errors_plot_between_cat(", "between_cat_list[['",i,"']]$sim_df_model_errors_sum", ")"
      )
    ), "\n\n"

    
    
    
  )
  
  qmd_lines <- c(qmd_lines, qmd_lines_i)
  
  
}

writeLines(qmd_lines, here::here("summary_qmds/between_cat.qmd"))