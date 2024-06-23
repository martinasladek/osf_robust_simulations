
# General ----------------------------------------

z <- function(x){(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)}

gen_df <- function(b0 = 0, b1 = 1, sd = 1, n = 1000){
  id = 1:n
  x1 = rnorm(n = n)
  e = rnorm(n = n, sd = sd)
  y = b0 + b1*x1
  data.frame(x1 = x1, y = y, e = e)
}

gen_df_sk_het <- function(n = 1000, 
                      b0 = 0, 
                      b1 = 0.25, 
                      b2 = 0.25, 
                      e) {
  
  id = 1:n
  x1 = rnorm(n = n)
  x2 = rnorm(n = n)
  
  e = e
  
  e1 = sample(e, size = length(e)/2, replace = FALSE)
  e2 = setdiff(e, e1)
  
  y = b0 + b1*x1 + b2*x2
  
  
  data.frame(x1 = x1, x2 = x2, y = y, e1 = e1, e2 = e2)
  
}

# Heteroscedasticity -----------------------------

vp1_val <- c(0)
vp2_val <- c(0.400, 0.600, 1.100, 1.400, 1.550)
vp3_val <- c(-0.200, -0.350, -0.650, -0.850, -0.950)
vp4_val <- c(0.200, 0.300, 0.500, 0.625, 0.700)


vp1 <- function(x1, vp_val) 1
vp2 <- function(x1, vp_val) (5 / (1+exp(-vp_val*(abs(x1)-1))))/2.5
vp3 <- function(x1, vp_val) exp(vp_val*abs(x1))
vp4 <- function(x1, vp_val) (5 / (1+exp(-vp_val*((x1)-1))))/2.5

# e.g.: 

modify_df <- function(df, vp, vp_val){
  df$y <- (df$y + vp(df$x1, vp_val)*df$e)
  df$x1 <- (df$x1)
  return(df)
}

modify_df_het_mult <- function(df, vp, vp_val){
  df$y <- df$y + (vp(df$x1, vp_val)*df$e1 + vp(df$x2, vp_val)*df$e2)
  df$x1 <- (df$x1)
  return(df)
}

# QLI -----------------------------

quantile_smoother	 <- function(y, x, 
                               prop_overlap = 0.75, # how much can the windows overlap
                               window_prop = 0.10, # what proportion of the sample size should each rolling window use?
                               tau = .95, # quantile 
                               window_alignment = c("center"), #
                               window_function = function(x) {quantile(x, tau)}
)
{
  
  sample_size <- length(y)
  window_size <- ceiling(sample_size*window_prop)
  
  window_distance <- window_size * (1-prop_overlap)
  
  
  # creating our new X and Y
  zoo.Y <- zoo(x = y, order.by = x)
  
  # center align 
  new.Y <- rollapply(zoo.Y,
                     width = window_size, 
                     FUN = window_function,
                     by = window_distance,
                     align = "center" 
  )
  
  new.X <- attributes(new.Y)$index
  new.Y <- as.numeric(new.Y) 
  
  
  # # lowess
  new.Y.mod <- lowess(new.Y~new.X)
  new.Y.loess <- new.Y.mod$y
  
  
  return(list(
    x = new.X, 
    y.loess = new.Y.loess
  ))
}

quantreg_interval_pred <- function(mod,
                                   predictor, 
                                   lower_quant = .025,
                                   upper_quant = .975,
                                   window_prop = 0.1,
                                   prop_overlap = 0.75
){
  
  q_lower <- quantile_smoother(
    y = mod$residuals,
    x = predictor,
    prop_overlap = prop_overlap,
    window_prop = window_prop,
    tau = lower_quant
  )
  
  q_upper <- quantile_smoother(
    y = mod$residuals,
    x = predictor,
    prop_overlap = prop_overlap,
    window_prop = window_prop,
    tau = upper_quant
  )
  
  quantreg_interval <- data.frame(
    x = q_lower$x,
    q_lower_y_loess = q_lower$y.loess,
    q_upper_y_loess = q_upper$y.loess,
    loess_wide = q_upper$y.loess - q_lower$y.loess
  )
  
  return(quantreg_interval)
}

quantreg_interval_pred_z <- function(mod,
                                   predictor, 
                                   lower_quant = .025,
                                   upper_quant = .975,
                                   window_prop = 0.1,
                                   prop_overlap = 0.75
){
  
  q_lower <- quantile_smoother(
    y = z(residuals(mod)),
    x = z(predictor),
    prop_overlap = prop_overlap,
    window_prop = window_prop,
    tau = lower_quant
  )
  
  q_upper <- quantile_smoother(
    y = z(residuals(mod)),
    x = z(predictor),
    prop_overlap = prop_overlap,
    window_prop = window_prop,
    tau = upper_quant
  )
  
  quantreg_interval <- data.frame(
    x = q_lower$x,
    q_lower_y_loess = q_lower$y.loess,
    q_upper_y_loess = q_upper$y.loess,
    loess_wide = q_upper$y.loess - q_lower$y.loess
  )
  
  return(quantreg_interval)
}

quantreg_interval_pred_z_e <- function(e,
                                       predictor, 
                                       lower_quant = .025,
                                       upper_quant = .975,
                                       window_prop = 0.1,
                                       prop_overlap = 0.75
){
  
  q_lower <- quantile_smoother(
    y = z(e),
    x = z(predictor),
    prop_overlap = prop_overlap,
    window_prop = window_prop,
    tau = lower_quant
  )
  
  q_upper <- quantile_smoother(
    y = z(e),
    x = z(predictor),
    prop_overlap = prop_overlap,
    window_prop = window_prop,
    tau = upper_quant
  )
  
  quantreg_interval <- data.frame(
    x = q_lower$x,
    q_lower_y_loess = q_lower$y.loess,
    q_upper_y_loess = q_upper$y.loess,
    loess_wide = q_upper$y.loess - q_lower$y.loess
  )
  
  return(quantreg_interval)
}

quantreg_mod_coefs <- function(quantreg_mod){
  quantreg_mod <- quantreg_mod$coefficients |> 
    t() |> 
    as.data.frame()
  
  names(quantreg_mod) <- paste0("quantreg_", names(quantreg_mod))
  
  return(quantreg_mod)
}

estimate_qli_x2_pos <- function(lm_mod, predictor){
  
  quantreg_df <- quantreg_interval_pred(mod = lm_mod, predictor = predictor)
  quantreg_mod <- lm(loess_wide ~ x + I(x^2) + I(x^3) + I(x^4), data = quantreg_df)
  quantreg_mod_coefs(quantreg_mod)[1, 3]
  
}

# Plotting ---------------------------------------

quick_dist <- function(x_arg, fill_arg = "darkcyan", colour_arg = "#005250", bins_arg = 60){
  ggplot2::ggplot(data = data.frame(), aes(x = x_arg)) + 
    geom_histogram(bins = bins_arg, fill = fill_arg, colour = colour_arg) + 
    theme_light()
}

quick_dens <- function(x_arg, fill_arg = "darkcyan", colour_arg = "#005250", bins_arg = 60){
  ggplot2::ggplot(data = data.frame(), aes(x = x_arg)) + 
    geom_density(fill = fill_arg, colour = colour_arg, alpha = 0.5, ) + 
    theme_light()
}

quick_scatter <- function(x, y){

    ggplot2::ggplot(data = data.frame(), aes(x = x, y = y)) + 
    geom_point(colour = "darkcyan", alpha = 0.1) + 
    theme_light()
  
}