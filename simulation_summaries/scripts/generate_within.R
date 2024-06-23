source(here::here("scripts/helpers.R"))


# set up files  -----------------------------------------------------------

within_cat <- all_files[grep("designs_within", all_files)]

within_cat_list <- purrr::map(
  .x = within_cat, 
  .f = readRDS
)

names(within_cat_list) <- sub(name_pattern, "", within_cat) |> 
  sub(x = _, replacement = "", pattern = ".rds") |> 
  toupper()


# processing tweaks -------------------------------------------------------

within_cat_list <- purrr::map(
  .x = within_cat_list, 
  .f = ~ purrr::map_at(
    .x, 
    .at = "sim_df_power_sum", 
    .f = ~ dplyr::mutate(
      .x, 
      name = factor(name, labels = c("OLS", "Bootstrap", "Trimming", "Trim + Bootstrap")))
  )
)


 
XW_list <- list(
  `CAT-2W` =  within_cat_list[["CAT-2W"]], 
  `CAT-3W` =  within_cat_list[["CAT-3W"]]
 )

XWXW_list <- setdiff(within_cat_list,  XW_list)

XW_list <- purrr::map(
  .x = XW_list,
  .f = ~ purrr::map_at(
    .x,
    .at = "sim_df_bias_sum",
    .f = ~ dplyr::mutate(
      .x,
      name = factor(name, labels = c("OLS", "Bootstrap", "Trimming", "Trim + Bootstrap"))
     
    )
  )
)

XW_list <- purrr::map(
  .x = XW_list,
  .f = ~ purrr::map_at(
    .x,
    .at = "sim_df_cover_sum",
    .f = ~ dplyr::mutate(
      .x,
      name = factor(name, labels = c("OLS", "Bootstrap", "Trimming", "Trim + Bootstrap"))
      
    )
  )
)

XW_list <- purrr::map(
  .x = XW_list,
  .f = ~ purrr::map_at(
    .x,
    .at = "sim_df_cover_shift_sum",
    .f = ~ dplyr::mutate(
      .x,
      name = factor(name, labels = c("OLS", "Bootstrap", "Trimming", "Trim + Bootstrap"))
      
    )
  )
)

XW_list <- purrr::map(
  .x = XW_list,
  .f = ~ purrr::map_at(
    .x,
    .at = "sim_df_ci_width_sum",
    .f = ~ dplyr::mutate(
      .x,
      name = factor(name, labels = c("OLS", "Bootstrap", "Trimming", "Trim + Bootstrap"))
      
    )
  )
)

XWXW_list <- purrr::map(
  .x = XWXW_list,
  .f = ~ purrr::map_at(
    .x,
    .at = "sim_df_bias_sum",
    .f = ~ dplyr::mutate(
      .x,
      name =factor(name, labels = c("OLS estimate", "Trimmed estimate"))
    )
  )
)

within_cat_list <- c(XW_list, XWXW_list) 

# generate qmd -----------------------------------------------------------


qmd_lines <- c()

for(i in names(within_cat_list)){
  
  design_name_i <- i
  
  unique_plot_groups_i <- within_cat_list[[i]]$sim_df_power_sum$plot_group |> unique() |> length()
  
  if(design_name_i %in% c("CAT-2W", "CAT-3W")){
    cat_23w_lines_i <- paste0(
      
      "### CI coverage", "\n\n",

      codechunk(
         settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
        code_string = paste0(
          "cover_plot_within_cat(", "within_cat_list[['",i,"']]$sim_df_cover_sum", ")"
        )
      ), "\n\n",
      
      "### CI bias-shifted coverage", "\n\n",

      codechunk(
         settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
        code_string = paste0(
          "cover_plot_within_cat(", "within_cat_list[['",i,"']]$sim_df_cover_shift_sum", ")"
        )
      ), "\n\n",
      
      "### CI width", "\n\n",

      codechunk(
        settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
        code_string = paste0(
          "ci_width_plot_within_cat(", "within_cat_list[['",i,"']]$sim_df_ci_width_sum", ")"
        )
      ), "\n\n"

      
    )
  } else {cat_23w_lines_i <- NULL}
  
  if(design_name_i == "CAT-3W"){
    cat_3w_lines_i <- paste0(
      
      "### Overall model fit", "\n\n",
      
      codechunk(
        settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
        code_string = paste0(
          "power_overall_plot_within_cat(", "within_cat_list[['",i,"']]$sim_df_overall_sum", ")"
        )
      ), "\n\n"
      
    )
  } else {
    cat_3w_lines_i <- NULL
  }
  
  qmd_lines_i <- paste0(
    
    "## ", design_name_i, "\n\n",
    
    cat_3w_lines_i,
    
    "### Power/false-positives", "\n\n",

    codechunk(
      settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
      code_string = paste0(
        "power_plot_within_cat(", "within_cat_list[['",i,"']]$sim_df_power_sum", ")"
      )
    ), "\n\n",
    
    "### Bias \n\n",
    
    "#### Mean bias \n\n",

    codechunk(
      settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
      code_string = paste0(
        "bias_mean_plot_within_cat(", "within_cat_list[['",i,"']]$sim_df_bias_sum", ")"
      )
    ), "\n\n",
    
    "#### SD bias \n\n",

    codechunk(
      settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
      code_string = paste0(
        "bias_sd_plot_within_cat(", "within_cat_list[['",i,"']]$sim_df_bias_sum", ")"
      )
    ), "\n\n",
    
    cat_23w_lines_i,
    
    "### Model errors", "\n\n",

    codechunk(
      settings_string = paste0("fig.height=", unique_plot_groups_i*1.3, ", fig.width=10"),
      code_string = paste0(
        "model_errors_plot_within_cat(", "within_cat_list[['",i,"']]$sim_df_model_errors_sum", ")"
      )
    ), "\n\n"
    
    
    
  )
  
  qmd_lines <- c(qmd_lines, qmd_lines_i)
  
}

writeLines(qmd_lines, here::here("summary_qmds/within_cat.qmd"))