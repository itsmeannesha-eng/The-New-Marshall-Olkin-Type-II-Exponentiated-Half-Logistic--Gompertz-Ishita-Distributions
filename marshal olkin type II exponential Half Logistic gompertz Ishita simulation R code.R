# Optimized & Fast Simulation Study for MO-T2EHL-Gom-Ishita Distribution
# Estimation Methods: Maximum Likelihood (ML), Least Squares (LS), WLS

# 1. Baseline CDF and PDF for Ishita Distribution
G_Ishita <- function(x, theta) {
  1 - (1 + (theta * x * (theta * x + 2)) / (theta^3 + 2)) * exp(-theta * x)
}

g_Ishita <- function(x, theta) {
  (theta^3 / (theta^3 + 2)) * (theta + x^2) * exp(-theta * x)
}

# 2. MO-T2EHL-Gom-Ishita Cumulative Distribution Function (CDF)
cdf_MO <- function(x, delta, alpha, gamma, theta) {
  G <- G_Ishita(x, theta)
  bar_G <- 1 - G
  
  inner_exp <- exp((1 / gamma) * (1 - bar_G^(-gamma)))
  brackets  <- inner_exp / (1 + (1 - inner_exp))
  
  delta_bar <- 1 - delta
  num <- delta * (brackets^alpha)
  den <- 1 - delta_bar * (brackets^alpha)
  
  return(1 - (num / den))
}

# 3. Fast Vectorized Random Variate Generation (Vectorized Bisection)
rmo_ishita_fast <- function(n, delta, alpha, gamma, theta, lower = 1e-7, upper = 50, tol_iter = 30) {
  u <- runif(n)
  delta_bar <- 1 - delta
  
  term1 <- ((1 - u) / (delta + (1 - u) * delta_bar))^(-1 / alpha)
  term2 <- 2 * ((term1 + 1)^(-1))
  term3 <- (1 - gamma * log(term2))^(-1 / gamma)
  
  target_G <- pmax(pmin(1 - term3, 1 - 1e-7), 1e-7)
  
  x_low <- rep(lower, n)
  x_high <- rep(upper, n)
  
  for (i in 1:tol_iter) {
    x_mid <- (x_low + x_high) / 2
    g_mid <- G_Ishita(x_mid, theta)
    
    left_side <- g_mid < target_G
    x_low[left_side] <- x_mid[left_side]
    x_high[!left_side] <- x_mid[!left_side]
  }
  
  return((x_low + x_high) / 2)
}

# 4. Objective Functions for Estimation Methods
neg_log_lik <- function(params, x) {
  delta <- params[1]; alpha <- params[2]; gamma <- params[3]; theta <- params[4]
  if (delta <= 0 || alpha <= 0 || gamma <= 0 || theta <= 0) return(1e10)
  
  n <- length(x)
  G <- G_Ishita(x, theta)
  g <- g_Ishita(x, theta)
  bar_G <- 1 - G
  
  if (any(bar_G <= 0) || any(G <= 0) || any(g <= 0)) return(1e10)
  
  inner_exp <- exp((1 / gamma) * (1 - bar_G^(-gamma)))
  ratio <- inner_exp / (1 + (1 - inner_exp))
  delta_bar <- 1 - delta
  
  log_lik <- n * log(2 * delta * alpha) +
    sum((alpha / gamma) * (1 - bar_G^(-gamma))) -
    (alpha + 1) * sum(log(1 + (1 - inner_exp))) -
    2 * sum(log(pmax(1e-10, 1 - delta_bar * (ratio^alpha)))) -
    (gamma + 1) * sum(log(bar_G)) +
    sum(log(g))
  
  if (is.na(log_lik) || is.infinite(log_lik)) return(1e10)
  return(-log_lik)
}

ls_obj <- function(params, x_sorted) {
  delta <- params[1]; alpha <- params[2]; gamma <- params[3]; theta <- params[4]
  if (delta <= 0 || alpha <= 0 || gamma <= 0 || theta <= 0) return(1e10)
  
  n <- length(x_sorted)
  i <- 1:n
  F_val <- cdf_MO(x_sorted, delta, alpha, gamma, theta)
  val <- sum((F_val - (i / (n + 1)))^2)
  if (is.na(val) || is.infinite(val)) return(1e10)
  return(val)
}

wls_obj <- function(params, x_sorted) {
  delta <- params[1]; alpha <- params[2]; gamma <- params[3]; theta <- params[4]
  if (delta <= 0 || alpha <= 0 || gamma <= 0 || theta <= 0) return(1e10)
  
  n <- length(x_sorted)
  i <- 1:n
  weights <- ((n + 1)^2 * (n + 2)) / (i * (n - i + 1))
  F_val <- cdf_MO(x_sorted, delta, alpha, gamma, theta)
  val <- sum(weights * (F_val - (i / (n + 1)))^2)
  if (is.na(val) || is.infinite(val)) return(1e10)
  return(val)
}

# 5. Main Simulation Routine
run_simulation <- function(true_params, sample_sizes, N = 500) {
  delta_true <- true_params["delta"]
  alpha_true <- true_params["alpha"]
  gamma_true <- true_params["gamma"]
  theta_true <- true_params["theta"]
  init_par   <- c(delta_true, alpha_true, gamma_true, theta_true)
  
  results_list <- list()
  
  for (n in sample_sizes) {
    cat(sprintf("Running simulation for n = %d ...\n", n))
    
    est_ML  <- matrix(NA, nrow = N, ncol = 4)
    est_LS  <- matrix(NA, nrow = N, ncol = 4)
    est_WLS <- matrix(NA, nrow = N, ncol = 4)
    
    set.seed(123)
    for (sim in 1:N) {
      x <- rmo_ishita_fast(n, delta_true, alpha_true, gamma_true, theta_true)
      x_sorted <- sort(x)
      
      # ML Estimation
      res_ML <- try(optim(init_par, neg_log_lik, x = x, method = "L-BFGS-B",
                          lower = c(1e-3, 1e-3, 1e-3, 1e-3)), silent = TRUE)
      if (!inherits(res_ML, "try-error") && res_ML$convergence == 0) {
        est_ML[sim, ] <- res_ML$par
      }
      
      # LS Estimation
      res_LS <- try(optim(init_par, ls_obj, x_sorted = x_sorted, method = "L-BFGS-B",
                          lower = c(1e-3, 1e-3, 1e-3, 1e-3)), silent = TRUE)
      if (!inherits(res_LS, "try-error") && res_LS$convergence == 0) {
        est_LS[sim, ] <- res_LS$par
      }
      
      # WLS Estimation
      res_WLS <- try(optim(init_par, wls_obj, x_sorted = x_sorted, method = "L-BFGS-B",
                           lower = c(1e-3, 1e-3, 1e-3, 1e-3)), silent = TRUE)
      if (!inherits(res_WLS, "try-error") && res_WLS$convergence == 0) {
        est_WLS[sim, ] <- res_WLS$par
      }
    }
    
    calc_metrics <- function(est_matrix) {
      metrics <- matrix(NA, nrow = 4, ncol = 2)
      colnames(metrics) <- c("ABIAS", "RMSE")
      rownames(metrics) <- c("delta", "alpha", "gamma", "theta")
      
      true_vec <- c(delta_true, alpha_true, gamma_true, theta_true)
      for (k in 1:4) {
        valid_est <- na.omit(est_matrix[, k])
        abias <- mean(valid_est) - true_vec[k]
        rmse  <- sqrt(mean((valid_est - true_vec[k])^2))
        metrics[k, ] <- c(abias, rmse)
      }
      return(metrics)
    }
    
    res_table <- data.frame(
      n = n,
      Parameter = c("delta", "alpha", "gamma", "theta"),
      ML_ABIAS  = calc_metrics(est_ML)[, "ABIAS"],
      ML_RMSE   = calc_metrics(est_ML)[, "RMSE"],
      LS_ABIAS  = calc_metrics(est_LS)[, "ABIAS"],
      LS_RMSE   = calc_metrics(est_LS)[, "RMSE"],
      WLS_ABIAS = calc_metrics(est_WLS)[, "ABIAS"],
      WLS_RMSE  = calc_metrics(est_WLS)[, "RMSE"]
    )
    
    results_list[[as.character(n)]] <- res_table
  }
  
  return(do.call(rbind, results_list))
}

# 6. Execute Simulation with New True Parameters
true_params <- c(delta = 0.8, alpha = 1.5, gamma = 0.5, theta = 2.0)
sample_sizes <- c(50, 100, 200, 500, 1000, 2000)

cat("Starting simulation with new parameter configuration...\n")
start_time <- Sys.time()

sim_results <- run_simulation(true_params, sample_sizes, N = 1000)

end_time <- Sys.time()
print(end_time - start_time)

# Print Final Summary Table
print(sim_results, digits = 4)
