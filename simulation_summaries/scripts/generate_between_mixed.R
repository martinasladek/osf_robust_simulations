source(here::here("scripts/helpers.R"))


# set up files  -----------------------------------------------------------

between_mixed <- all_files[grep("designs_mixed_predictors", all_files)] |> 
  grep(x = _, pattern = "split|preds", invert = TRUE, value = TRUE)

between_mixed_list <- purrr::map(
  .x = between_mixed, 
  .f = readRDS
)

names(between_mixed_list) <- sub(name_pattern, "", between_mixed) |> 
  sub(x = _, replacement = "", pattern = ".rds")

# 
# names(between_mixed_list$`N-NUM-CAT-INT_split`) <- names(between_mixed_list$`N-NUM-CAT-INT_split`) |> 
#   sub(x = _, replacement = "", pattern = "sim_export_N-NUM-CAT-INT_preds_") |> 
#   sub(x = _, replacement = " predictors", pattern = ".rds")
# 
# 
# names(between_mixed_list$`N-NUM-CAT_split`) <- names(between_mixed_list$`N-NUM-CAT_split`) |> 
#   sub(x = _, replacement = "", pattern = "sim_export_N-NUM-CAT_preds_") |> 
#   sub(x = _, replacement = " predictors", pattern = ".rds")
# processing tweaks -------------------------------------------------------

# add nice names to power plots 

between_mixed_list <- purrr::map(
  .x = between_mixed_list, 
  .f = ~ purrr::map_at(
    .x, 
    .at = "sim_df_power_sum", 
    .f = ~ dplyr::mutate(
      .x, 
      name = factor(name, labels = c("OLS", "Bootstrap", "HC4", "MM/KS", "MM")))
  )
)


# add nice names to power plots 

between_mixed_list <- purrr::map(
  .x = between_mixed_list, 
  .f = ~ purrr::map_at(
    .x, 
    .at = "sim_df_cover_sum", 
    .f = ~ dplyr::mutate(
      .x, 
      name = factor(name, labels = c("OLS", "Bootstrap", "HC4", "MM/KS", "MM")))
  )
)

# add nice names to error plots 

between_mixed_list <- purrr::map(
  .x = between_mixed_list, 
  .f = ~ purrr::map_at(
    .x, 
    .at = "sim_df_model_errors_sum", 
    .f = ~ dplyr::mutate(
      .x, 
      name = factor(name, labels = c("OLS", "Bootstrap", "HC4", "MM/KS", "MM")))
  )
)

# tidy up overall power name and plot group 

between_mixed_list <- purrr::map(
  .x = between_mixed_list, 
  .f = ~ purrr::map_at(
    .x, 
    .at = "sim_df_overall_sum", 
    .f = ~ dplyr::mutate(
      .x, 
      name = "OLS"
    )
  )
)

between_mixed_list <- purrr::map(
  .x = between_mixed_list, 
  .f = ~ purrr::map_at(
    .x, 
    .at = "sim_df_overall_sum", 
    .f = ~ dplyr::mutate(
      .x, 
      plot_group = stringr::str_replace_all(string = plot_group, pattern = "\nKurtosis", "; Kurtosis"), 
      het = stringr::str_replace_all(string = het, pattern = "; het_shape", "\nhet_shape"),
      plot_group = paste0(plot_group, "\n", het)
    )
  )
)

# generate qmd -----------------------------------------------------------


qmd_lines <- c()

#for(i in names(between_mixed_list)){

for(i in names(between_mixed_list)){
  for(j in c(1,3,8)){
    
    design_name_i <- i
    
    unique_plot_groups_i <-( between_mixed_list[[i]]$sim_df_overall_sum$plot_group |> unique() |> length())/3
    
    qmd_lines_i <- paste0(
      
      "## ", design_name_i, ": ", j, " continuous predictor(s)", "\n\n",
      
      "### Overall model fit", "\n\n",
      
      codechunk(
        settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
        code_string = paste0(
          "power_overall_plot_between_cat(", "between_mixed_list[['",i,"']]$sim_df_overall_sum |> dplyr::filter(n_cont_preds ==",j,")", ")"
        )
      ), "\n\n",
      
      "### Power/false-positives", "\n\n",
      
      codechunk(
        settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
        code_string = paste0(
          "power_plot_between_cat(", "between_mixed_list[['",i,"']]$sim_df_power_sum |> dplyr::filter(n_cont_preds ==",j,")", ")"
        )
      ), "\n\n",
      
      
      "### Bias \n\n",
      
      "#### Mean bias \n\n",
      
      codechunk(
        settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
        code_string = paste0(
          "bias_mean_plot_between_cat(", "between_mixed_list[['",i,"']]$sim_df_bias_sum |> dplyr::filter(n_cont_preds ==",j,")", ")"
        )
      ), "\n\n",
      
      "#### SD bias \n\n",
      
      codechunk(
        settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
        code_string = paste0(
          "bias_sd_plot_between_cat(", "between_mixed_list[['",i,"']]$sim_df_bias_sum |> dplyr::filter(n_cont_preds ==",j,")", ")"
        )
      ), "\n\n",
      
      "### CI coverage", "\n\n",
      
      codechunk(
        settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
        code_string = paste0(
          "cover_plot_between_cat(", "between_mixed_list[['",i,"']]$sim_df_cover_sum |> dplyr::filter(n_cont_preds ==",j,")", ")"
        )
      ), "\n\n",
      
      "### CI bias-shifted coverage", "\n\n",
      
      codechunk(
        settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
        code_string = paste0(
          "cover_plot_between_cat(", "between_mixed_list[['",i,"']]$sim_df_cover_shift_sum |> dplyr::filter(n_cont_preds ==",j,")", ")"
        )
      ), "\n\n",
      
      "### CI width", "\n\n",
      
      codechunk(
        settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
        code_string = paste0(
          "ci_width_plot_between_cat(", "between_mixed_list[['",i,"']]$sim_df_ci_width_sum |> dplyr::filter(n_cont_preds ==",j,")", ")"
        )
      ), "\n\n",
      
      
      "### Model errors", "\n\n",
      
      codechunk(
        settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
        code_string = paste0(
          "model_errors_plot_between_cat(", "between_mixed_list[['",i,"']]$sim_df_model_errors_sum |> dplyr::filter(n_cont_preds ==",j,")", ")"
        )
      ), "\n\n"
      
      
      
      
    )
    
    qmd_lines <- c(qmd_lines, qmd_lines_i)
    
    
  }
  
  writeLines(
    qmd_lines, 
    here::here(
      paste0("summary_qmds/between_mixed_",design_name_i, ".qmd")
    )
  )
  
  qmd_lines <- c()
  
}
