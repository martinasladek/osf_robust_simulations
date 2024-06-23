

# file set up -------------------------------------------------------------

name_pattern <- "^.*sim_export_"

all_files <- list.files("../../", pattern = "sim_export*", recursive = TRUE, full.names = TRUE)



# file.copy(
#   from = all_files, 
#   to = "/Users/ms2027/Library/Mobile Documents/com~apple~CloudDocs/Documents/uni/uni_phd/writing/robust_methods_thesis/chapter_6_simulations/r_docs/manuscripts/data/processed_data"
# )

# general helpers ---------------------------------------------------------

codechunk <- function(code_string, settings_string = ""){
  paste0(
    "```{r, ", settings_string, "}", "\n", 
    code_string, "\n", 
    "```", "\n\n"
  )
}

# between_cat plotting helpers --------------------------------------------

power_overall_plot_between_cat <- function(overall_power_sum){
  
  power_plot <- overall_power_sum |> 
    dplyr::mutate(name_helper = paste0(name, " [", prop_sig ,"]")) |> 
    ggplot2::ggplot(aes(
      x = n_bw, y = prop_sig, linetype = name, colour = name, shape = name, fill = name)) + 
    geom_line() + 
    ggiraph::geom_point_interactive(size = 2, aes(tooltip = name_helper)) + 
    scale_x_continuous(breaks = unique(overall_power_sum$n_bw)) + 
    labs(x = "Sample size", y = "Proportion of p < .05") + 
    facet_grid(
      rows = vars(plot_group), cols = vars(b)) + 
    scale_colour_viridis_d(option = "A", end = 0.8) + 
    scale_fill_viridis_d(option = "A", end = 0.8) + 
    scale_linetype_manual(values = c("solid", "dotted", "dashed", "dotdash", "twodash")) + 
    scale_shape_manual(values = c(21,22,23,24,25)) + 
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.text.y.right = element_text(angle = 0), 
      panel.spacing = unit(1, "lines")
    ) 
  
  ggiraph::girafe(ggobj = power_plot)
  
  
}


power_plot_between_cat <- function(power_sum){
  
  power_plot <- power_sum |> 
    dplyr::mutate(name_helper = paste0(name, " [", prop_sig ,"]")) |> 
    ggplot2::ggplot(aes(
      x = n_bw, y = prop_sig, linetype = name, colour = name, shape = name, fill = name)) + 
    geom_line() + 
    ggiraph::geom_point_interactive(size = 2, aes(tooltip = name_helper)) + 
    scale_x_continuous(breaks = unique(power_sum$n_bw)) + 
    labs(x = "Sample size", y = "Proportion of p < .05") + 
    facet_grid(
      rows = vars(plot_group), cols = vars(b)) + 
    scale_colour_viridis_d(option = "A", end = 0.8) + 
    scale_fill_viridis_d(option = "A", end = 0.8) + 
    scale_linetype_manual(values = c("solid", "dotted", "dashed", "dotdash", "twodash")) + 
    scale_shape_manual(values = c(21,22,23,24,25)) + 
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.text.y.right = element_text(angle = 0), 
      panel.spacing = unit(1, "lines")
    ) 
  
  ggiraph::girafe(ggobj = power_plot)
  
  
}



bias_mean_plot_between_cat <- function(bias_sum){
  
  bias_mean_plot <- bias_sum |> 
    dplyr::mutate(name_helper = paste0(name, " [", mean_diff ,"]")) |> 
    ggplot2::ggplot(aes(x = n_bw, y = mean_diff, 
                        linetype = name, colour = name, shape = name, fill = name)) +
    geom_line() + 
    ggiraph::geom_point_interactive(size = 2, aes(tooltip = name_helper)) + 
    scale_x_continuous(breaks = unique(bias_sum$n_bw)) + 
    labs(x = "Sample size", y = "Mean difference") + 
    facet_grid(
      rows = vars(plot_group), cols = vars(b),
      scales = "free_y") + 
    scale_colour_viridis_d(option = "A", end = 0.8) + 
    scale_fill_viridis_d(option = "A", end = 0.8) + 
    scale_linetype_manual(values = c("solid", "dotted", "dotdash", "twodash")) + 
    scale_shape_manual(values = c(21,22,24,25)) + 
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.text.y.right = element_text(angle = 0), 
      panel.spacing = unit(1, "lines")
    ) 
  
  ggiraph::girafe(ggobj = bias_mean_plot)
  
  
}

bias_sd_plot_between_cat <- function(bias_sum){
  
  bias_sd_plot <- bias_sum |> 
    dplyr::mutate(name_helper = paste0(name, " [", sd_diff ,"]")) |> 
    ggplot2::ggplot(aes(x = n_bw, y = sd_diff, 
                        linetype = name, colour = name, shape = name, fill = name)) +
    geom_line() + 
    ggiraph::geom_point_interactive(size = 2, aes(tooltip = name_helper)) + 
    scale_x_continuous(breaks = unique(bias_sum$n_bw)) + 
    labs(x = "Sample size", y = "SD") + 
    facet_grid(
      rows = vars(plot_group), cols = vars(b),
      scales = "free_y") + 
    scale_colour_viridis_d(option = "A", end = 0.8) + 
    scale_fill_viridis_d(option = "A", end = 0.8) + 
    scale_linetype_manual(values = c("solid", "dotted", "dotdash", "twodash")) + 
    scale_shape_manual(values = c(21,22,24,25)) + 
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.text.y.right = element_text(angle = 0), 
      panel.spacing = unit(1, "lines")
    )
  
  ggiraph::girafe(ggobj = bias_sd_plot)
  
  
}


cover_plot_between_cat <- function(cover_sum){
  
  cover_plot <- cover_sum |> 
    dplyr::mutate(name_helper = paste0(name, " [", prop_covered ,"]")) |> 
    ggplot2::ggplot(aes(x = n_bw, y = prop_covered, 
                        linetype = name, colour = name, shape = name, fill = name)) + 
    geom_line() + 
    ggiraph::geom_point_interactive(size = 2, aes(tooltip = name_helper)) + 
    scale_x_continuous(breaks = unique(cover_sum$n_bw)) + 
    labs(x = "Sample size", y = "Coverage") + 
    facet_grid(
      rows = vars(plot_group), cols = vars(b),
      scales = "free_y") + 
    scale_colour_viridis_d(option = "A", end = 0.8) + 
    scale_fill_viridis_d(option = "A", end = 0.8) + 
    scale_linetype_manual(values = c("solid", "dotted", "dashed", "dotdash", "twodash")) + 
    scale_shape_manual(values = c(21,22,23,24,25)) + 
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.text.y.right = element_text(angle = 0), 
      panel.spacing = unit(1, "lines")
    ) 
  
  ggiraph::girafe(ggobj = cover_plot)
  
  
}

ci_width_plot_between_cat <- function(ci_width_sum){
  
  ci_width_plot <- ci_width_sum |> 
    dplyr::mutate(name_helper = paste0(name, " [", mean_width ,"]")) |> 
    ggplot2::ggplot(aes(x = n_bw, y = mean_width, 
                        linetype = name, colour = name, shape = name, fill = name)) + 
    geom_line() + 
    ggiraph::geom_point_interactive(size = 2, aes(tooltip = name_helper)) + 
    scale_x_continuous(breaks = unique(ci_width_sum$n_bw)) + 
    labs(x = "Sample size", y = "CI width") + 
    facet_grid(
      rows = vars(plot_group), cols = vars(b),
      scales = "free_y") + 
    scale_colour_viridis_d(option = "A", end = 0.8) + 
    scale_fill_viridis_d(option = "A", end = 0.8) + 
    scale_linetype_manual(values = c("solid", "dotted", "dashed", "dotdash", "twodash")) + 
    scale_shape_manual(values = c(21,22,23,24,25)) + 
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.text.y.right = element_text(angle = 0), 
      panel.spacing = unit(1, "lines")
    ) 
  
  ggiraph::girafe(ggobj = ci_width_plot)
  
  
}


model_errors_plot_between_cat <- function(model_errors_sum){
  
  model_errors_plot <- model_errors_sum |> 
    dplyr::mutate(name_helper = paste0(name, " [", n_fail ,"]")) |> 
    ggplot2::ggplot(aes(x = n_bw, y = n_fail, 
                        linetype = name, colour = name, shape = name, fill = name))+
    geom_line() + 
    ggiraph::geom_point_interactive(size = 2, aes(tooltip = name_helper)) + 
    scale_x_continuous(breaks = unique(model_errors_sum$n_bw)) + 
    labs(x = "Sample size", y = "Proportion of failed models") + 
    facet_grid(
      rows = vars(plot_group), cols = vars(b),
      scales = "free_y") + 
    scale_colour_viridis_d(option = "A", end = 0.8) + 
    scale_fill_viridis_d(option = "A", end = 0.8) + 
    scale_linetype_manual(values = c("solid", "dotted", "dashed", "dotdash", "twodash")) + 
    scale_shape_manual(values = c(21,22,23,24,25)) + 
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.text.y.right = element_text(angle = 0), 
      panel.spacing = unit(1, "lines")
    ) 
  
  ggiraph::girafe(ggobj = model_errors_plot)
  
  
}
  


# within_cat plotting helpers ---------------------------------------------

power_overall_plot_within_cat <- function(overall_power_sum){
  
  power_plot <- overall_power_sum |> 
    dplyr::mutate(
      name = factor(
        name, 
        levels = c("afex_nc_p", "afex_gg_p"), 
        labels = c("OLS (uncorrected)", "OLS (Greenhouse-Geisser)")
      ),
      name_helper = paste0(name, " [", prop_sig ,"]")) |> 
    ggplot2::ggplot(aes(
      x = n_rm, y = prop_sig, linetype = name, colour = name, shape = name, fill = name)) + 
    geom_line() + 
    ggiraph::geom_point_interactive(size = 2, aes(tooltip = name_helper)) + 
    scale_x_continuous(breaks = unique(overall_power_sum$n_rm)) + 
    labs(x = "Sample size", y = "Proportion of p < .05") + 
    facet_grid(
      rows = vars(plot_group), cols = vars(b)) + 
    scale_colour_viridis_d(option = "A", end = 0.8) + 
    scale_fill_viridis_d(option = "A", end = 0.8) + 
    scale_linetype_manual(values = c("solid", "dotted", "dashed", "dotdash", "twodash")) + 
    scale_shape_manual(values = c(21,22,23,24,25)) + 
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.text.y.right = element_text(angle = 0), 
      panel.spacing = unit(1, "lines")
    ) 
  
  ggiraph::girafe(ggobj = power_plot)
  
  
}


power_plot_within_cat <- function(power_sum){
  
  
  power_plot <- power_sum |> 
    dplyr::mutate(name_helper = paste0(name, " [", prop_sig ,"]")) |> 
    ggplot2::ggplot(
      aes(
        x = n_rm, y = prop_sig, linetype = name, colour = name, shape = name, fill = name 
      )
    )+ 
    geom_line() + 
    ggiraph::geom_point_interactive(size = 2, aes(tooltip = name_helper)) + 
    scale_x_continuous(breaks = unique(power_sum$n_rm)) + 
    labs(x = "Sample size", y = "Proportion of p < .05") + 
    facet_grid(
      rows = vars(plot_group), cols = vars(b)
    ) + 
    scale_colour_viridis_d(option = "A", end = 0.8) + 
    scale_fill_viridis_d(option = "A", end = 0.8) + 
    scale_linetype_manual(values = c("solid", "dotted", "dashed", "dotdash", "twodash")) + 
    scale_shape_manual(values = c(21,22,23,24,25)) + 
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.text.y.right = element_text(angle = 0), 
      panel.spacing = unit(1, "lines")
    ) 
  
  ggiraph::girafe(ggobj = power_plot)
  
  

  
  
}


bias_mean_plot_within_cat <- function(bias_sum){
  
  bias_mean_plot <- bias_sum |> 
    dplyr::mutate(name_helper = paste0(name, " [", mean_diff ,"]")) |> 
    ggplot2::ggplot(aes(x = n_rm, y = mean_diff, 
                        linetype = name, colour = name, shape = name, fill = name)) +
    geom_line() + 
    ggiraph::geom_point_interactive(size = 2, aes(tooltip = name_helper)) + 
    scale_x_continuous(breaks = unique(bias_sum$n_rm)) + 
    labs(x = "Sample size", y = "Mean difference") + 
    facet_grid(
      rows = vars(plot_group), cols = vars(b),
      scales = "free_y") + 
    scale_colour_viridis_d(option = "A", end = 0.8) + 
    scale_fill_viridis_d(option = "A", end = 0.8) + 
    scale_linetype_manual(values = c("solid", "dotted", "dotdash", "twodash")) + 
    scale_shape_manual(values = c(21,22,24,25)) + 
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.text.y.right = element_text(angle = 0), 
      panel.spacing = unit(1, "lines")
    ) 
  
  ggiraph::girafe(ggobj = bias_mean_plot)
  
  
}

bias_sd_plot_within_cat <- function(bias_sum){
  
  bias_sd_plot <- bias_sum |> 
    dplyr::mutate(name_helper = paste0(name, " [", sd_diff ,"]")) |> 
    ggplot2::ggplot(aes(x = n_rm, y = sd_diff, 
                        linetype = name, colour = name, shape = name, fill = name)) +
    geom_line() + 
    ggiraph::geom_point_interactive(size = 2, aes(tooltip = name_helper)) + 
    scale_x_continuous(breaks = unique(bias_sum$n_rm)) + 
    labs(x = "Sample size", y = "SD") + 
    facet_grid(
      rows = vars(plot_group), cols = vars(b),
      scales = "free_y") + 
    scale_colour_viridis_d(option = "A", end = 0.8) + 
    scale_fill_viridis_d(option = "A", end = 0.8) + 
    scale_linetype_manual(values = c("solid", "dotted", "dotdash", "twodash")) + 
    scale_shape_manual(values = c(21,22,24,25)) + 
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.text.y.right = element_text(angle = 0), 
      panel.spacing = unit(1, "lines")
    )
  
  ggiraph::girafe(ggobj = bias_sd_plot)
  
  
}

model_errors_plot_within_cat <- function(model_errors_sum){
  
  model_errors_plot <- model_errors_sum |> 
    dplyr::mutate(name_helper = paste0(name, " [", n_fail ,"]")) |> 
    ggplot2::ggplot(aes(x = n_rm, y = n_fail, 
                        linetype = name, colour = name, shape = name, fill = name))+
    geom_line() + 
    ggiraph::geom_point_interactive(size = 2, aes(tooltip = name_helper)) + 
    scale_x_continuous(breaks = unique(model_errors_sum$n_rm)) + 
    labs(x = "Sample size", y = "Proportion of failed models") + 
    facet_grid(
      rows = vars(plot_group), cols = vars(b),
      scales = "free_y") + 
    scale_colour_viridis_d(option = "A", end = 0.8) + 
    scale_fill_viridis_d(option = "A", end = 0.8) + 
    scale_linetype_manual(values = c("solid", "dotted", "dashed", "dotdash", "twodash")) + 
    scale_shape_manual(values = c(21,22,23,24,25)) + 
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.text.y.right = element_text(angle = 0), 
      panel.spacing = unit(1, "lines")
    ) 
  
  ggiraph::girafe(ggobj = model_errors_plot)
  
  
}


cover_plot_within_cat <- function(cover_sum){
  
  cover_plot <- cover_sum |> 
    dplyr::mutate(name_helper = paste0(name, " [", prop_covered ,"]")) |> 
    ggplot2::ggplot(aes(x = n_rm, y = prop_covered, 
                        linetype = name, colour = name, shape = name, fill = name)) + 
    geom_line() + 
    ggiraph::geom_point_interactive(size = 2, aes(tooltip = name_helper)) + 
    scale_x_continuous(breaks = unique(cover_sum$n_rm)) + 
    labs(x = "Sample size", y = "Coverage") + 
    facet_grid(
      rows = vars(plot_group), cols = vars(b),
      scales = "free_y") + 
    scale_colour_viridis_d(option = "A", end = 0.8) + 
    scale_fill_viridis_d(option = "A", end = 0.8) + 
    scale_linetype_manual(values = c("solid", "dotted", "dashed", "dotdash", "twodash")) + 
    scale_shape_manual(values = c(21,22,23,24,25)) + 
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.text.y.right = element_text(angle = 0), 
      panel.spacing = unit(1, "lines")
    ) 
  
  ggiraph::girafe(ggobj = cover_plot)
  
  
}

ci_width_plot_within_cat <- function(ci_width_sum){
  
  ci_width_plot <- ci_width_sum |> 
    dplyr::mutate(name_helper = paste0(name, " [", mean_width ,"]")) |> 
    ggplot2::ggplot(aes(x = n_rm, y = mean_width, 
                        linetype = name, colour = name, shape = name, fill = name)) + 
    geom_line() + 
    ggiraph::geom_point_interactive(size = 2, aes(tooltip = name_helper)) + 
    scale_x_continuous(breaks = unique(ci_width_sum$n_rm)) + 
    labs(x = "Sample size", y = "CI width") + 
    facet_grid(
      rows = vars(plot_group), cols = vars(b),
      scales = "free_y") + 
    scale_colour_viridis_d(option = "A", end = 0.8) + 
    scale_fill_viridis_d(option = "A", end = 0.8) + 
    scale_linetype_manual(values = c("solid", "dotted", "dashed", "dotdash", "twodash")) + 
    scale_shape_manual(values = c(21,22,23,24,25)) + 
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.text.y.right = element_text(angle = 0), 
      panel.spacing = unit(1, "lines")
    ) 
  
  ggiraph::girafe(ggobj = ci_width_plot)
  
  
}