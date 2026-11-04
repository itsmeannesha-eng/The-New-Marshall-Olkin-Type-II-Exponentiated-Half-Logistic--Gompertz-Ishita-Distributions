# R Script: Bayesian Computation, Tables, and Advanced Plots for the
# MO-T2EHL-Gom-Ishita Distribution

# Install/Load required packages
required_packages <- c("coda", "ggplot2", "gridExtra", "MASS", "KernSmooth")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 1. Baseline and Model Probability Functions

G_ishita <- function(x, theta) {
  1 - ((exp(-theta * x)) / (theta^3 + 2)) * (theta^3 + 2 + theta^2 * x + theta * x^2 + x^3)
}

g_ishita <- function(x, theta) {
  (theta^3 / (theta^3 + 2)) * (theta + x^2) * exp(-theta * x)
}

pdf_mo_t2ehl_gom_ishita <- function(x, delta, alpha, gamma, theta) {
  G <- G_ishita(x, theta)
  g_val <- g_ishita(x, theta)
  bar_G <- 1 - G
  bar_delta <- 1 - delta
  
  # Avoid underflow/overflow bounds
  bar_G[bar_G <= 1e-15] <- 1e-15
  
  A_val <- exp((1 / gamma) * (1 - bar_G^(-gamma)))
  term <- A_val / (1 + (1 - A_val))
  term[term <= 0] <- 1e-15
  term[term >= 1] <- 1 - 1e-15
  
  f_val <- 2 * delta * alpha * exp((alpha / gamma) * (1 - bar_G^(-gamma))) *
    bar_G^(-gamma - 1) * (1 + (1 - A_val))^(-alpha - 1) *
    g_val * (1 - bar_delta * term^alpha)^(-2)
  
  return(f_val)
}

log_likelihood <- function(params, x_data) {
  delta <- params[1]; alpha <- params[2]; gamma <- params[3]; theta <- params[4]
  if (delta <= 0 || alpha <= 0 || gamma <= 0 || theta <= 0) return(-Inf)
  
  densities <- pdf_mo_t2ehl_gom_ishita(x_data, delta, alpha, gamma, theta)
  if (any(densities <= 0) || any(!is.finite(densities))) return(-Inf)
  
  sum(log(densities))
}


# 2. Prior Specification and Log-Posterior

log_prior <- function(params, hyper_params) {
  a <- hyper_params$a
  b <- hyper_params$b
  
  lp_delta <- dgamma(params[1], shape = a[1], rate = b[1], log = TRUE)
  lp_alpha <- dgamma(params[2], shape = a[2], rate = b[2], log = TRUE)
  lp_gamma <- dgamma(params[3], shape = a[3], rate = b[3], log = TRUE)
  lp_theta <- dgamma(params[4], shape = a[4], rate = b[4], log = TRUE)
  
  if (!all(is.finite(c(lp_delta, lp_alpha, lp_gamma, lp_theta)))) return(-Inf)
  return(lp_delta + lp_alpha + lp_gamma + lp_theta)
}

log_posterior <- function(params, x_data, hyper_params) {
  ll <- log_likelihood(params, x_data)
  if (ll == -Inf) return(-Inf)
  lp <- log_prior(params, hyper_params)
  if (lp == -Inf) return(-Inf)
  return(ll + lp)
}


# 3. Metropolis-Hastings (MH) MCMC Sampler

run_mcmc_mh <- function(x_data, init_params, sigma_mle, M = 12000, M_0 = 2000, 
                        hyper_params = list(a = c(0.001, 0.001, 0.001, 0.001), 
                                            b = c(0.001, 0.001, 0.001, 0.001))) {
  c_scale <- 2.38 / sqrt(4)
  prop_cov <- (c_scale^2) * sigma_mle
  
  samples <- matrix(NA, nrow = M, ncol = 4)
  colnames(samples) <- c("delta", "alpha", "gamma", "theta")
  
  current_params <- init_params
  current_phi <- log(current_params)
  current_log_post <- log_posterior(current_params, x_data, hyper_params)
  
  accepted_count <- 0
  
  for (t in 1:M) {
    prop_phi <- MASS::mvrnorm(n = 1, mu = current_phi, Sigma = prop_cov)
    prop_params <- exp(prop_phi)
    
    prop_log_post <- log_posterior(prop_params, x_data, hyper_params)
    
    if (is.finite(prop_log_post)) {
      # Jacobian adjustment for log-normal proposal transformation
      log_alpha_mh <- (prop_log_post + sum(prop_phi)) - (current_log_post + sum(current_phi))
      
      if (log(runif(1)) < log_alpha_mh) {
        current_params <- prop_params
        current_phi <- prop_phi
        current_log_post <- prop_log_post
        if (t > M_0) accepted_count <- accepted_count + 1
      }
    }
    samples[t, ] <- current_params
  }
  
  cat("MCMC Sampling Completed. Acceptance Rate:", round(accepted_count / (M - M_0) * 100, 2), "%\n")
  return(samples[(M_0 + 1):M, ])
}


# 4. Simulation Execution (Using Synthetic Data Example)

set.seed(2026)
# Simulate sample observations from a baseline distribution profile
x_sample <- rgamma(60, shape = 2, rate = 1.2)

# Initial MLE starting points and covariance approximation matrix
init_mle <- c(delta = 0.5, alpha = 1.2, gamma = 0.8, theta = 1.5)
mle_cov_matrix <- diag(c(0.015, 0.04, 0.025, 0.03))

# Run the MCMC Chain
posterior_chain <- run_mcmc_mh(x_data = x_sample, init_params = init_mle, 
                               sigma_mle = mle_cov_matrix, M = 15000, M_0 = 3000)


# 5. Summary Table Generation (Point & Interval Estimates)

generate_bayesian_table <- function(post_samples) {
  post_mean <- colMeans(post_samples)
  post_median <- apply(post_samples, 2, median)
  hpd_intervals <- coda::HPDinterval(coda::as.mcmc(post_samples), prob = 0.95)
  ess_vals <- apply(post_samples, 2, coda::effectiveSize)
  
  summary_table <- data.frame(
    Parameter = c("delta", "alpha", "gamma", "theta"),
    Posterior_Mean = round(post_mean, 4),
    Posterior_Median = round(post_median, 4),
    HPD_Lower_95 = round(hpd_intervals[, "lower"], 4),
    HPD_Upper_95 = round(hpd_intervals[, "upper"], 4),
    Effective_Sample_Size = round(ess_vals, 1)
  )
  rownames(summary_table) <- NULL
  return(summary_table)
}

# Re-run table generation
bayes_table <- generate_bayesian_table(posterior_chain)
print("=== TABLE: Bayesian Estimation Summary Results ===")
print(bayes_table)


# 6. Colorful & Advanced Diagnostic Plots (ggplot2)

# Convert matrix to data frame for ggplot visualization
df_mcmc <- data.frame(posterior_chain)
df_mcmc$Iteration <- 1:nrow(df_mcmc)

# A. Trace Plots
p_delta_trace <- ggplot(df_mcmc, aes(x = Iteration, y = delta)) +
  geom_line(color = "#2E86C1", alpha = 0.8) +
  labs(title = "Trace Plot: delta", x = "Iteration", y = expression(delta)) +
  theme_minimal(base_size = 12)

p_alpha_trace <- ggplot(df_mcmc, aes(x = Iteration, y = alpha)) +
  geom_line(color = "#28B463", alpha = 0.8) +
  labs(title = "Trace Plot: alpha", x = "Iteration", y = expression(alpha)) +
  theme_minimal(base_size = 12)

p_gamma_trace <- ggplot(df_mcmc, aes(x = Iteration, y = gamma)) +
  geom_line(color = "#D4AC0D", alpha = 0.8) +
  labs(title = "Trace Plot: gamma", x = "Iteration", y = expression(gamma)) +
  theme_minimal(base_size = 12)

p_theta_trace <- ggplot(df_mcmc, aes(x = Iteration, y = theta)) +
  geom_line(color = "#CB4335", alpha = 0.8) +
  labs(title = "Trace Plot: theta", x = "Iteration", y = expression(theta)) +
  theme_minimal(base_size = 12)

# Arrange Trace Plots into a Grid
grid.arrange(p_delta_trace, p_alpha_trace, p_gamma_trace, p_theta_trace, ncol = 2,
             top = "MCMC Trace Plots for Model Parameters")

# B. Marginal Posterior Density Plots with HPD Shading Intervals
plot_marginal_density <- function(data_vector, param_name, color_hex) {
  hpd <- coda::HPDinterval(coda::as.mcmc(data_vector), prob = 0.95)
  df_dens <- data.frame(val = data_vector)
  
  ggplot(df_dens, aes(x = val)) +
    geom_density(fill = color_hex, alpha = 0.4, color = color_hex, linewidth = 1) +
    geom_vline(xintercept = mean(data_vector), linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_vline(xintercept = hpd[1], linetype = "dotted", color = "darkblue", linewidth = 0.8) +
    geom_vline(xintercept = hpd[2], linetype = "dotted", color = "darkblue", linewidth = 0.8) +
    labs(title = paste("Posterior Density:", param_name), x = param_name, y = "Density") +
    theme_minimal(base_size = 12)
}

p_delta_dens <- plot_marginal_density(df_mcmc$delta, "delta", "#2E86C1")
p_alpha_dens <- plot_marginal_density(df_mcmc$alpha, "alpha", "#28B463")
p_gamma_dens <- plot_marginal_density(df_mcmc$gamma, "gamma", "#D4AC0D")
p_theta_dens <- plot_marginal_density(df_mcmc$theta, "theta", "#CB4335")

grid.arrange(p_delta_dens, p_alpha_dens, p_gamma_dens, p_theta_dens, ncol = 2,
             top = "Marginal Posterior Densities with 95% HPD Bounds")

# C. Advanced Bivariate Contour Plots (delta vs alpha AND gamma vs theta)

# 1. Contour Plot: delta vs alpha
p_contour_1 <- ggplot(df_mcmc, aes(x = delta, y = alpha)) +
  geom_point(color = "#8E44AD", alpha = 0.1, size = 0.8) +
  geom_density_2d(color = "#1B4F72", linewidth = 0.8) +
  labs(title = "Bivariate Contour: delta vs alpha",
       x = expression(delta), y = expression(alpha)) +
  theme_minimal(base_size = 12)

# 2. Contour Plot: gamma vs theta
p_contour_2 <- ggplot(df_mcmc, aes(x = gamma, y = theta)) +
  geom_point(color = "#E67E22", alpha = 0.1, size = 0.8) +
  geom_density_2d(color = "#78281F", linewidth = 0.8) +
  labs(title = "Bivariate Contour: gamma vs theta",
       x = expression(gamma), y = expression(theta)) +
  theme_minimal(base_size = 12)

# Combine both plots into a side-by-side grid layout
grid.arrange(p_contour_1, p_contour_2, ncol = 2,
             top = "Bivariate Posterior Densities and Contour Maps")
