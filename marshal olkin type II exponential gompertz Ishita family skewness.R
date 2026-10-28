if (!require("plotly")) install.packages("plotly")
library(plotly)


# Baseline Ishita & Distribution Functions

g_ishita <- function(x, theta) {
  ifelse(x <= 0, 0, (theta^3 / (theta^3 + 2)) * (theta + x^2) * exp(-theta * x))
}

G_ishita <- function(x, theta) {
  ifelse(x <= 0, 0, 1 - (1 + (theta * x * (theta * x + 2)) / (theta^3 + 2)) * exp(-theta * x))
}

pdf_MOT2EHL_Gom_Ishita <- function(x, delta, alpha, gamma, theta) {
  res <- numeric(length(x))
  pos <- x > 0
  if (!any(pos)) return(res)
  
  x_pos <- x[pos]
  d_bar <- 1 - delta
  
  G_val <- G_ishita(x_pos, theta)
  g_val <- g_ishita(x_pos, theta)
  S_bar <- 1 - G_val
  
  # Protect against numerical overflow/underflow in exponents
  S_bar <- pmax(S_bar, 1e-12)
  
  A_x <- 1 - (S_bar)^(-gamma)
  B_x <- exp((1 / gamma) * A_x)
  C_x <- B_x / (2 - B_x)
  
  term1 <- 2 * delta * alpha
  term2 <- exp((alpha / gamma) * A_x)
  term3 <- S_bar^(-gamma - 1)
  term4 <- (2 - B_x)^(-(alpha + 1))
  term5 <- g_val
  term6 <- (1 - d_bar * (C_x)^alpha)^(-2)
  
  val <- term1 * term2 * term3 * term4 * term5 * term6
  res[pos] <- ifelse(is.finite(val), val, 0)
  return(res)
}

#  Skewness 
compute_skewness <- function(delta, alpha, gamma, theta) {
  # upper bound for integration to avoid infinite tail failure
  upper_bound <- 50 
  
  raw_moment <- function(k) {
    res <- tryCatch({
      integrate(function(x) (x^k) * pdf_MOT2EHL_Gom_Ishita(x, delta, alpha, gamma, theta),
                lower = 1e-7, upper = upper_bound, subdivisions = 500, rel.tol = 1e-4)$value
    }, error = function(e) NA_real_)
    return(res)
  }
  
  mu1 <- raw_moment(1)
  mu2 <- raw_moment(2)
  mu3 <- raw_moment(3)
  
  if (is.na(mu1) || is.na(mu2) || is.na(mu3)) return(NA_real_)
  
  var_x <- mu2 - mu1^2
  if (is.na(var_x) || var_x <= 0) return(NA_real_)
  
  sigma <- sqrt(var_x)
  mu3_central <- mu3 - 3 * mu1 * mu2 + 2 * mu1^3
  
  skew <- mu3_central / (sigma^3)
  return(if (is.finite(skew)) skew else NA_real_)
}


# Grid Evaluation
gamma_fixed <- 1.0
theta_fixed <- 1.2

alpha_seq <- seq(0.5, 2.5, length.out = 15)
delta_seq <- seq(0.5, 2.5, length.out = 15)

z_skewness <- matrix(NA_real_, nrow = length(alpha_seq), ncol = length(delta_seq))

cat("Computing surface values...\n")
for (i in seq_along(alpha_seq)) {
  for (j in seq_along(delta_seq)) {
    z_skewness[i, j] <- compute_skewness(delta = delta_seq[j], 
                                         alpha = alpha_seq[i], 
                                         gamma = gamma_fixed, 
                                         theta = theta_fixed)
  }
}

z_matrix <- t(z_skewness)

if (all(is.na(z_matrix))) {
  stop("Integration failed for all parameter values. Adjust alpha/delta ranges.")
} else if (any(is.na(z_matrix))) {
  z_matrix[is.na(z_matrix)] <- mean(z_matrix, na.rm = TRUE)
}

#Plot
fig <- plot_ly(x = alpha_seq, y = delta_seq, z = z_matrix, type = "surface")
fig <- plotly::layout(
  fig,
  title = paste0(" "),
  scene = list(
    xaxis = list(title = "alpha "),
    yaxis = list(title = "delta "),
    zaxis = list(title = "Skewness")
  )
)

fig


# Install the required R package
if (!requireNamespace("reticulate", quietly = TRUE)) {
  install.packages("reticulate")
}

# Export the plot as a PNG image
plotly::save_image(
  fig,
  file = "Skewness_3D_Plot.png",
  width = 1600,
  height = 1200,
  scale = 2
)





if (!require("plotly")) install.packages("plotly")
library(plotly)


# Ishita & Distribution Functions
g_ishita <- function(x, theta) {
  ifelse(x <= 0, 0, (theta^3 / (theta^3 + 2)) * (theta + x^2) * exp(-theta * x))
}

G_ishita <- function(x, theta) {
  ifelse(x <= 0, 0, 1 - (1 + (theta * x * (theta * x + 2)) / (theta^3 + 2)) * exp(-theta * x))
}

pdf_MOT2EHL_Gom_Ishita <- function(x, delta, alpha, gamma, theta) {
  res <- numeric(length(x))
  pos <- x > 0
  if (!any(pos)) return(res)
  
  x_pos <- x[pos]
  d_bar <- 1 - delta
  
  G_val <- G_ishita(x_pos, theta)
  g_val <- g_ishita(x_pos, theta)
  S_bar <- 1 - G_val
  
  # Protect against numerical overflow/underflow in exponents
  S_bar <- pmax(S_bar, 1e-12)
  
  A_x <- 1 - (S_bar)^(-gamma)
  B_x <- exp((1 / gamma) * A_x)
  C_x <- B_x / (2 - B_x)
  
  term1 <- 2 * delta * alpha
  term2 <- exp((alpha / gamma) * A_x)
  term3 <- S_bar^(-gamma - 1)
  term4 <- (2 - B_x)^(-(alpha + 1))
  term5 <- g_val
  term6 <- (1 - d_bar * (C_x)^alpha)^(-2)
  
  val <- term1 * term2 * term3 * term4 * term5 * term6
  res[pos] <- ifelse(is.finite(val), val, 0)
  return(res)
}

# Skewness 
compute_skewness <- function(delta, alpha, gamma, theta) {
  # Practical upper bound for integration to avoid infinite tail failure
  upper_bound <- 50 
  
  raw_moment <- function(k) {
    res <- tryCatch({
      integrate(function(x) (x^k) * pdf_MOT2EHL_Gom_Ishita(x, delta, alpha, gamma, theta),
                lower = 1e-7, upper = upper_bound, subdivisions = 500, rel.tol = 1e-4)$value
    }, error = function(e) NA_real_)
    return(res)
  }
  
  mu1 <- raw_moment(1)
  mu2 <- raw_moment(2)
  mu3 <- raw_moment(3)
  
  if (is.na(mu1) || is.na(mu2) || is.na(mu3)) return(NA_real_)
  
  var_x <- mu2 - mu1^2
  if (is.na(var_x) || var_x <= 0) return(NA_real_)
  
  sigma <- sqrt(var_x)
  mu3_central <- mu3 - 3 * mu1 * mu2 + 2 * mu1^3
  
  skew <- mu3_central / (sigma^3)
  return(if (is.finite(skew)) skew else NA_real_)
}

# 3. Grid Evaluation
# Fixed Parameters 
delta_fixed <- 1.0
alpha_fixed <- 1.5

# Variable Ranges (x-axis: theta, y-axis: gamma)
theta_seq <- seq(0.5, 2.5, length.out = 15)
gamma_seq <- seq(0.5, 2.5, length.out = 15)

z_skewness <- matrix(NA_real_, nrow = length(theta_seq), ncol = length(gamma_seq))

cat("Computing surface values...\n")
for (i in seq_along(theta_seq)) {
  for (j in seq_along(gamma_seq)) {
    z_skewness[i, j] <- compute_skewness(delta = delta_fixed, 
                                         alpha = alpha_fixed, 
                                         gamma = gamma_seq[j], 
                                         theta = theta_seq[i])
  }
}

z_matrix <- t(z_skewness)

# Fill any residual NA values via surface mean to protect color scale
if (all(is.na(z_matrix))) {
  stop("Integration failed for all parameter values. Adjust theta/gamma ranges.")
} else if (any(is.na(z_matrix))) {
  z_matrix[is.na(z_matrix)] <- mean(z_matrix, na.rm = TRUE)
}

# 4. Interactive Plotly 3D Surface Plot
fig <- plot_ly(x = theta_seq, y = gamma_seq, z = z_matrix, type = "surface")
fig <- plotly::layout(
  fig,
  title = paste0(" "),
  scene = list(
    xaxis = list(title = "theta (θ)"),
    yaxis = list(title = "gamma (γ)"),
    zaxis = list(title = "Skewness")
  )
)

fig

# Confirm the saved file location
cat(
  "Plot saved at:",
  normalizePath(
    file.path(save_folder, "Skewness_Theta_Gamma_3D.html"),
    winslash = "/",
    mustWork = FALSE
  ),
  "\n"
)


