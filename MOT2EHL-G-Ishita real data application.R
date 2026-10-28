# 1. MOTOR DATA

library(boot)

data("motor", package = "boot")

x <- na.omit(boot::motor$times)
x <- as.numeric(x[is.finite(x) & x > 0])
n <- length(x)

cat("Sample size =", n, "\n")
print(summary(x))

# 2. BASELINE ISHITA PDF AND CDF

dIshita <- function(x, theta) {
  ans <- numeric(length(x))
  ok <- is.finite(x) & x > 0

  ans[ok] <- theta^3 / (theta^3 + 2) *
    (theta + x[ok]^2) * exp(-theta * x[ok])

  ans
}

pIshita <- function(x, theta) {
  ans <- numeric(length(x))
  ans[x <= 0] <- 0
  ans[is.infinite(x) & x > 0] <- 1

  ok <- is.finite(x) & x > 0
  z <- theta * x[ok]

  ans[ok] <- theta^3 / (theta^3 + 2) *
    ((1 - exp(-z)) +
       2 / theta^3 *
       (1 - exp(-z) * (1 + z + z^2 / 2)))

  pmin(1, pmax(0, ans))
}


# 3. ISHITA-POLYNOMIAL PDF AND CDF
# f_m(x) = C_m(theta) *
#          (theta + x^2 + sum_{k=3}^m x^k) exp(-theta*x)

Cm <- function(theta, m) {
  terms <- if (m >= 3) {
    sapply(3:m, function(k) factorial(k) / theta^(k + 1))
  } else {
    numeric(0)
  }

  1 / (1 + 2 / theta^3 + sum(terms))
}

dPolyIshita <- function(x, theta, m) {
  ans <- numeric(length(x))
  ok <- is.finite(x) & x > 0

  if (any(ok)) {
    xx <- x[ok]
    poly <- theta + xx^2

    if (m >= 3) {
      for (k in 3:m) poly <- poly + xx^k
    }

    ans[ok] <- Cm(theta, m) * poly * exp(-theta * xx)
  }

  ans
}

pPolyIshita <- function(x, theta, m) {
  ans <- numeric(length(x))
  ans[x <= 0] <- 0
  ans[is.infinite(x) & x > 0] <- 1

  ok <- is.finite(x) & x > 0

  if (any(ok)) {
    z <- theta * x[ok]

    # Integral of theta * exp(-theta*t)
    val <- 1 - exp(-z)

    # Integral of t^2 * exp(-theta*t)
    val <- val + 2 / theta^3 * pgamma(z, shape = 3)

    # Integrals of t^k * exp(-theta*t), k = 3,...,m
    if (m >= 3) {
      for (k in 3:m) {
        val <- val +
          factorial(k) / theta^(k + 1) *
          pgamma(z, shape = k + 1)
      }
    }

    ans[ok] <- Cm(theta, m) * val
  }

  pmin(1, pmax(0, ans))
}


# 4. PROPOSED MOT2EHL-GOM-ISHITA MODEL
# Parameters: delta, alpha, gamma, theta

pProposed <- function(x, delta, alpha, gamma, theta) {
  G <- pIshita(x, theta)
  S <- pmax(1 - G, 1e-14)

  v <- (1 - S^(-gamma)) / gamma
  q <- exp(pmax(v, -700))
  A <- (q / (2 - q))^alpha

  F <- 1 - delta * A / (1 - (1 - delta) * A)

  F[x <= 0] <- 0
  F[is.infinite(x) & x > 0] <- 1

  pmin(1, pmax(0, F))
}

dProposed <- function(x, delta, alpha, gamma, theta) {
  G <- pIshita(x, theta)
  g <- dIshita(x, theta)
  S <- pmax(1 - G, 1e-14)

  v <- (1 - S^(-gamma)) / gamma
  q <- exp(pmax(v, -700))
  A <- (q / (2 - q))^alpha

  den <- 1 - (1 - delta) * A

  f <- delta * 2 * alpha * A *
    S^(-gamma - 1) * g /
    ((2 - q) * den^2)

  f[!is.finite(f) | x <= 0] <- 0
  pmax(f, 0)
}


# 5. MAXIMUM-LIKELIHOOD ESTIMATION

# Fit Ishita model
nll_ishita <- function(par) {
  theta <- exp(par[1])
  d <- dIshita(x, theta)

  if (any(!is.finite(d)) || any(d <= 0)) return(1e100)

  -sum(log(d))
}

fitI <- optim(
  log(1 / mean(x)),
  nll_ishita,
  method = "BFGS",
  hessian = TRUE
)

thetaI <- exp(fitI$par)
llI <- -fitI$value

# Fit polynomial Ishita models for m = 3, 4, 5
fitPoly <- function(m) {
  obj <- function(par) {
    theta <- exp(par[1])
    d <- dPolyIshita(x, theta, m)

    if (any(!is.finite(d)) || any(d <= 0)) return(1e100)

    -sum(log(d))
  }

  fit <- optim(
    log(1 / mean(x)),
    obj,
    method = "BFGS"
  )

  list(
    theta = exp(fit$par),
    logLik = -fit$value,
    convergence = fit$convergence
  )
}

polyFits <- lapply(3:5, fitPoly)
names(polyFits) <- paste0("Ishita-Polynomial-m", 3:5)

# Fit proposed model using log-transformed positive parameters
nll_proposed <- function(par) {
  pars <- exp(par)

  delta <- pars[1]
  alpha <- pars[2]
  gamma <- pars[3]
  theta <- pars[4]

  d <- dProposed(x, delta, alpha, gamma, theta)

  if (any(!is.finite(d)) || any(d <= 0)) return(1e100)

  -sum(log(d))
}

fitP <- optim(
  par = log(c(1, 1, 1, 1 / mean(x))),
  fn = nll_proposed,
  method = "Nelder-Mead",
  control = list(maxit = 50000, reltol = 1e-9)
)

proposedPars <- exp(fitP$par)
llP <- -fitP$value

cat("\nIshita theta estimate:\n")
print(thetaI)

cat("\nPolynomial estimates:\n")
print(lapply(polyFits, function(z) z$theta))

cat("\nProposed model estimates (delta, alpha, gamma, theta):\n")
print(proposedPars)

cat("\nProposed model convergence code:", fitP$convergence, "\n")


# 6. MODEL COMPARISON: LOG-LIKELIHOOD, AIC, BIC, KS

modelNames <- c(
  "Ishita",
  names(polyFits),
  "MOT2EHL-Gom-Ishita"
)

logLikValues <- c(
  llI,
  sapply(polyFits, function(z) z$logLik),
  llP
)

numPars <- c(1, 1, 1, 1, 4)

# Fitted CDFs for KS statistic
cdfFunctions <- c(
  list(function(z) pIshita(z, thetaI)),
  lapply(3:5, function(m) {
    local({
      mm <- m
      tt <- polyFits[[paste0("Ishita-Polynomial-m", mm)]]$theta
      function(z) pPolyIshita(z, tt, mm)
    })
  }),
  list(function(z) pProposed(
    z,
    proposedPars[1],
    proposedPars[2],
    proposedPars[3],
    proposedPars[4]
  ))
)

ksStatistic <- function(data, cdfFun) {
  z <- sort(data)
  nn <- length(z)
  Fz <- cdfFun(z)

  max(
    max(seq_len(nn) / nn - Fz),
    max(Fz - (seq_len(nn) - 1) / nn)
  )
}

KS <- sapply(cdfFunctions, function(ff) ksStatistic(x, ff))

# Approximate KS p-values for each fitted CDF
# Note: ordinary KS p-values are approximate because model
# parameters have been estimated from the same sample.
KS_p_value <- sapply(cdfFunctions, function(ff) {
  tryCatch(
    suppressWarnings(ks.test(x, ff)$p.value),
    error = function(e) NA_real_
  )
})

comparison <- data.frame(
  Model = modelNames,
  Parameters = numPars,
  LogLik = logLikValues,
  AIC = -2 * logLikValues + 2 * numPars,
  BIC = -2 * logLikValues + log(n) * numPars,
  KS = KS,
  KS_p_value = KS_p_value
)

comparison <- comparison[order(comparison$AIC), ]
rownames(comparison) <- NULL

cat("\nMODEL COMPARISON TABLE\n")
print(comparison, digits = 5)

write.csv(
  comparison,
  "Ishita_model_comparison_motor.csv",
  row.names = FALSE
)


# 7. HISTOGRAM WITH FITTED DENSITIES
# ADJUSTED X-AXIS AND Y-AXIS LIMITS

# Extended x-grid for smooth fitted density curves
xgrid <- seq(
  0,
  max(x) * 1.15,
  length.out = 2000
)

# Calculate fitted density values
dens_I <- dIshita(xgrid, thetaI)

dens_P3 <- dPolyIshita(
  xgrid,
  polyFits[["Ishita-Polynomial-m3"]]$theta,
  3
)

dens_P4 <- dPolyIshita(
  xgrid,
  polyFits[["Ishita-Polynomial-m4"]]$theta,
  4
)

dens_P5 <- dPolyIshita(
  xgrid,
  polyFits[["Ishita-Polynomial-m5"]]$theta,
  5
)

dens_Prop <- dProposed(
  xgrid,
  proposedPars[1],
  proposedPars[2],
  proposedPars[3],
  proposedPars[4]
)

# Calculate histogram heights without plotting
hist_data <- hist(
  x,
  breaks = "FD",
  plot = FALSE
)

# Find the largest finite density across histogram and models
all_heights <- c(
  hist_data$density,
  dens_I,
  dens_P3,
  dens_P4,
  dens_P5,
  dens_Prop
)

all_heights <- all_heights[is.finite(all_heights)]

ymax <- max(all_heights, na.rm = TRUE)

# Add 15% space above the highest value
ylim_max <- ymax * 1.15

# Plot histogram with adjusted axis limits
hist(
  x,
  probability = TRUE,
  breaks = "FD",
  col = "grey90",
  border = "white",
  xlim = c(0, max(x) * 1.15),
  ylim = c(0, ylim_max),
  main = "",
  xlab = "Motor testing time",
  ylab = "Density",
  cex.main = 1.3,
  cex.lab = 1.2
)

# Add all fitted density curves using solid colored lines
lines(
  xgrid, dens_I,
  col = "blue",
  lwd = 2.5,
  lty = 1
)

lines(
  xgrid, dens_P3,
  col = "red",
  lwd = 2.5,
  lty = 1
)

lines(
  xgrid, dens_P4,
  col = "darkgreen",
  lwd = 2.5,
  lty = 1
)

lines(
  xgrid, dens_P5,
  col = "purple",
  lwd = 2.5,
  lty = 1
)

lines(
  xgrid, dens_Prop,
  col = "orange",
  lwd = 2.5,
  lty = 1
)

# Legend
legend(
  "topright",
  legend = c(
    "Ishita",
    "Polynomial m=3",
    "Polynomial m=4",
    "Polynomial m=5",
    "Proposed MOT2EHL-Gom-Ishita"
  ),
  col = c(
    "blue",
    "red",
    "darkgreen",
    "purple",
    "orange"
  ),
  lty = 1,
  lwd = 2.5,
  cex = 0.8,
  bty = "n",
  inset = 0.02
)


# 8. EMPIRICAL CDF AND FITTED CDF OVERLAYS

plot(
  ecdf(x),
  verticals = TRUE,
  do.points = FALSE,
  main = "",
  xlab = "Motor testing time",
  ylab = "Cumulative probability",
  xlim = c(0, max(x) * 1.15),
  ylim = c(0, 1)
)

lines(
  xgrid,
  pIshita(xgrid, thetaI),
  col = "blue",
  lwd = 2.5,
  lty = 1
)

for (m in 3:5) {
  tt <- polyFits[[paste0("Ishita-Polynomial-m", m)]]$theta

  lines(
    xgrid,
    pPolyIshita(xgrid, tt, m),
    col = c("red", "darkgreen", "purple")[m - 2],
    lwd = 2.5,
    lty = 1
  )
}

lines(
  xgrid,
  pProposed(
    xgrid,
    proposedPars[1],
    proposedPars[2],
    proposedPars[3],
    proposedPars[4]
  ),
  col = "orange",
  lwd = 2.5,
  lty = 1
)

legend(
  "bottomright",
  legend = c(
    "Empirical CDF",
    "Ishita",
    "Polynomial m=3",
    "Polynomial m=4",
    "Polynomial m=5",
    "Proposed MOT2EHL-Gom-Ishita"
  ),
  col = c(
    "black",
    "blue",
    "red",
    "darkgreen",
    "purple",
    "orange"
  ),
  lty = 1,
  lwd = c(2.5, 2.5, 2.5, 2.5, 2.5, 2.5),
  cex = 0.8,
  bty = "n"
)


# 9. OUTPUT FILE LOCATION
cat("\nOutput directory:\n")
print(getwd())

cat("\nCSV comparison table saved as:\n")
cat("Ishita_model_comparison_motor.csv\n")










# 2. GEYSER DATA

if (!requireNamespace("MASS", quietly = TRUE)) {
  install.packages("MASS")
}

library(MASS)

data("geyser", package = "MASS")

x <- MASS::geyser$waiting
x <- as.numeric(x[is.finite(x) & x > 0])
n <- length(x)

cat("Sample size =", n, "\n")
print(summary(x))


# 2. BASELINE ISHITA PDF AND CDF

dIshita <- function(x, theta) {
  ans <- numeric(length(x))
  ok <- is.finite(x) & x > 0
  
  ans[ok] <- theta^3 / (theta^3 + 2) *
    (theta + x[ok]^2) * exp(-theta * x[ok])
  
  ans
}

pIshita <- function(x, theta) {
  ans <- numeric(length(x))
  ans[x <= 0] <- 0
  ans[is.infinite(x) & x > 0] <- 1
  
  ok <- is.finite(x) & x > 0
  z <- theta * x[ok]
  
  ans[ok] <- theta^3 / (theta^3 + 2) *
    ((1 - exp(-z)) +
       2 / theta^3 *
       (1 - exp(-z) * (1 + z + z^2 / 2)))
  
  pmin(1, pmax(0, ans))
}


# 3. ISHITA-POLYNOMIAL PDF AND CDF
# f_m(x) = C_m(theta) *
#          (theta + x^2 + sum_{k=3}^m x^k) exp(-theta*x)

Cm <- function(theta, m) {
  terms <- if (m >= 3) {
    sapply(3:m, function(k) factorial(k) / theta^(k + 1))
  } else {
    numeric(0)
  }
  
  1 / (1 + 2 / theta^3 + sum(terms))
}

dPolyIshita <- function(x, theta, m) {
  ans <- numeric(length(x))
  ok <- is.finite(x) & x > 0
  
  if (any(ok)) {
    xx <- x[ok]
    poly <- theta + xx^2
    
    if (m >= 3) {
      for (k in 3:m) {
        poly <- poly + xx^k
      }
    }
    
    ans[ok] <- Cm(theta, m) * poly * exp(-theta * xx)
  }
  
  ans
}

pPolyIshita <- function(x, theta, m) {
  ans <- numeric(length(x))
  ans[x <= 0] <- 0
  ans[is.infinite(x) & x > 0] <- 1
  
  ok <- is.finite(x) & x > 0
  
  if (any(ok)) {
    z <- theta * x[ok]
    
    # Integral of theta * exp(-theta*t)
    val <- 1 - exp(-z)
    
    # Integral of t^2 * exp(-theta*t)
    val <- val + 2 / theta^3 * pgamma(z, shape = 3)
    
    # Integrals of t^k * exp(-theta*t), k = 3,...,m
    if (m >= 3) {
      for (k in 3:m) {
        val <- val +
          factorial(k) / theta^(k + 1) *
          pgamma(z, shape = k + 1)
      }
    }
    
    ans[ok] <- Cm(theta, m) * val
  }
  
  pmin(1, pmax(0, ans))
}


# 4. PROPOSED MOT2EHL-GOM-ISHITA MODEL
# Parameters: delta, alpha, gamma, theta

pProposed <- function(x, delta, alpha, gamma, theta) {
  G <- pIshita(x, theta)
  S <- pmax(1 - G, 1e-14)
  
  v <- (1 - S^(-gamma)) / gamma
  q <- exp(pmax(v, -700))
  A <- (q / (2 - q))^alpha
  
  F <- 1 - delta * A / (1 - (1 - delta) * A)
  
  F[x <= 0] <- 0
  F[is.infinite(x) & x > 0] <- 1
  
  pmin(1, pmax(0, F))
}

dProposed <- function(x, delta, alpha, gamma, theta) {
  G <- pIshita(x, theta)
  g <- dIshita(x, theta)
  S <- pmax(1 - G, 1e-14)
  
  v <- (1 - S^(-gamma)) / gamma
  q <- exp(pmax(v, -700))
  A <- (q / (2 - q))^alpha
  
  den <- 1 - (1 - delta) * A
  
  f <- delta * 2 * alpha * A *
    S^(-gamma - 1) * g /
    ((2 - q) * den^2)
  
  f[!is.finite(f) | x <= 0] <- 0
  pmax(f, 0)
}


# 5. MAXIMUM-LIKELIHOOD ESTIMATION

# Fit baseline Ishita model
nll_ishita <- function(par) {
  theta <- exp(par[1])
  d <- dIshita(x, theta)
  
  if (any(!is.finite(d)) || any(d <= 0)) {
    return(1e100)
  }
  
  -sum(log(d))
}

fitI <- optim(
  log(1 / mean(x)),
  nll_ishita,
  method = "BFGS",
  hessian = TRUE
)

thetaI <- exp(fitI$par)
llI <- -fitI$value


# Fit polynomial Ishita models for m = 3, 4, 5
fitPoly <- function(m) {
  obj <- function(par) {
    theta <- exp(par[1])
    d <- dPolyIshita(x, theta, m)
    
    if (any(!is.finite(d)) || any(d <= 0)) {
      return(1e100)
    }
    
    -sum(log(d))
  }
  
  fit <- optim(
    log(1 / mean(x)),
    obj,
    method = "BFGS"
  )
  
  list(
    theta = exp(fit$par),
    logLik = -fit$value,
    convergence = fit$convergence
  )
}

polyFits <- lapply(3:5, fitPoly)
names(polyFits) <- paste0("Ishita-Polynomial-m", 3:5)


# Fit proposed model using log-transformed positive parameters
nll_proposed <- function(par) {
  pars <- exp(par)
  
  delta <- pars[1]
  alpha <- pars[2]
  gamma <- pars[3]
  theta <- pars[4]
  
  d <- dProposed(x, delta, alpha, gamma, theta)
  
  if (any(!is.finite(d)) || any(d <= 0)) {
    return(1e100)
  }
  
  -sum(log(d))
}

fitP <- optim(
  par = log(c(1, 1, 1, 1 / mean(x))),
  fn = nll_proposed,
  method = "Nelder-Mead",
  control = list(maxit = 50000, reltol = 1e-9)
)

proposedPars <- exp(fitP$par)
llP <- -fitP$value


# Display parameter estimates
cat("\nIshita theta estimate:\n")
print(thetaI)

cat("\nPolynomial estimates:\n")
print(lapply(polyFits, function(z) z$theta))

cat("\nProposed model estimates (delta, alpha, gamma, theta):\n")
print(proposedPars)

cat("\nProposed model convergence code:", fitP$convergence, "\n")



# 6. MODEL COMPARISON: LOGLIKELIHOOD, AIC, BIC, KS AND P-VALUE

modelNames <- c(
  "Ishita",
  names(polyFits),
  "MOT2EHL-Gom-Ishita"
)

logLikValues <- c(
  llI,
  sapply(polyFits, function(z) z$logLik),
  llP
)

numPars <- c(1, 1, 1, 1, 4)


# Fitted CDFs for KS statistic and p-value
cdfFunctions <- c(
  list(function(z) pIshita(z, thetaI)),
  
  lapply(3:5, function(m) {
    local({
      mm <- m
      tt <- polyFits[[paste0("Ishita-Polynomial-m", mm)]]$theta
      function(z) pPolyIshita(z, tt, mm)
    })
  }),
  
  list(function(z) {
    pProposed(
      z,
      proposedPars[1],
      proposedPars[2],
      proposedPars[3],
      proposedPars[4]
    )
  })
)


# Calculate KS statistic
ksStatistic <- function(data, cdfFun) {
  z <- sort(data)
  nn <- length(z)
  Fz <- cdfFun(z)
  
  max(
    max(seq_len(nn) / nn - Fz),
    max(Fz - (seq_len(nn) - 1) / nn)
  )
}

KS <- sapply(
  cdfFunctions,
  function(ff) ksStatistic(x, ff)
)


# Approximate KS p-values
# These are approximate because model parameters are estimated
# using the same sample. Parametric bootstrap is preferable
# for formal goodness-of-fit testing.
KS_p_value <- sapply(cdfFunctions, function(ff) {
  tryCatch(
    suppressWarnings(ks.test(x, ff)$p.value),
    error = function(e) NA_real_
  )
})


# Create model comparison table
comparison <- data.frame(
  Model = modelNames,
  Parameters = numPars,
  LogLik = logLikValues,
  AIC = -2 * logLikValues + 2 * numPars,
  BIC = -2 * logLikValues + log(n) * numPars,
  KS = KS,
  KS_p_value = KS_p_value
)

comparison <- comparison[order(comparison$AIC), ]
rownames(comparison) <- NULL

cat("\nMODEL COMPARISON TABLE\n")
print(comparison, digits = 5)

write.csv(
  comparison,
  "Ishita_model_comparison_geyser.csv",
  row.names = FALSE
)


# 7. HISTOGRAM WITH FITTED DENSITIES
# ADJUSTED X-AXIS AND Y-AXIS LIMITS

# Extended x-grid for smooth fitted density curves
xgrid <- seq(
  0,
  max(x) * 1.15,
  length.out = 2000
)


# Calculate fitted density values
dens_I <- dIshita(xgrid, thetaI)

dens_P3 <- dPolyIshita(
  xgrid,
  polyFits[["Ishita-Polynomial-m3"]]$theta,
  3
)

dens_P4 <- dPolyIshita(
  xgrid,
  polyFits[["Ishita-Polynomial-m4"]]$theta,
  4
)

dens_P5 <- dPolyIshita(
  xgrid,
  polyFits[["Ishita-Polynomial-m5"]]$theta,
  5
)

dens_Prop <- dProposed(
  xgrid,
  proposedPars[1],
  proposedPars[2],
  proposedPars[3],
  proposedPars[4]
)


# Calculate histogram heights without plotting
hist_data <- hist(
  x,
  breaks = "FD",
  plot = FALSE
)


# Find the largest finite density across histogram and models
all_heights <- c(
  hist_data$density,
  dens_I,
  dens_P3,
  dens_P4,
  dens_P5,
  dens_Prop
)

all_heights <- all_heights[is.finite(all_heights)]

ymax <- max(all_heights, na.rm = TRUE)
ylim_max <- ymax * 1.15


# Plot histogram with adjusted axis limits
hist(
  x,
  probability = TRUE,
  breaks = "FD",
  col = "grey90",
  border = "white",
  xlim = c(0, max(x) * 1.15),
  ylim = c(0, ylim_max),
  main = "",
  xlab = "Waiting time (minutes)",
  ylab = "Density",
  cex.main = 1.3,
  cex.lab = 1.2
)


# Add fitted density curves
lines(
  xgrid, dens_I,
  col = "blue",
  lwd = 2.5,
  lty = 1
)

lines(
  xgrid, dens_P3,
  col = "red",
  lwd = 2.5,
  lty = 1
)

lines(
  xgrid, dens_P4,
  col = "darkgreen",
  lwd = 2.5,
  lty = 1
)

lines(
  xgrid, dens_P5,
  col = "purple",
  lwd = 2.5,
  lty = 1
)

lines(
  xgrid, dens_Prop,
  col = "orange",
  lwd = 2.5,
  lty = 1
)


# Legend
legend(
  "topleft",
  legend = c(
    "Ishita",
    "Polynomial m=3",
    "Polynomial m=4",
    "Polynomial m=5",
    "Proposed MOT2EHL-Gom-Ishita"
  ),
  col = c(
    "blue",
    "red",
    "darkgreen",
    "purple",
    "orange"
  ),
  lty = 1,
  lwd = 2.5,
  cex = 0.8,
  bty = "n",
  inset = 0.02
)


# 8. EMPIRICAL CDF AND FITTED CDF OVERLAYS

plot(
  ecdf(x),
  verticals = TRUE,
  do.points = FALSE,
  main = "",
  xlab = "Waiting time (minutes)",
  ylab = "Cumulative probability",
  xlim = c(0, max(x) * 1.15),
  ylim = c(0, 1)
)


# Ishita CDF
lines(
  xgrid,
  pIshita(xgrid, thetaI),
  col = "blue",
  lwd = 2.5,
  lty = 1
)


# Polynomial Ishita CDFs
for (m in 3:5) {
  tt <- polyFits[[paste0("Ishita-Polynomial-m", m)]]$theta
  
  lines(
    xgrid,
    pPolyIshita(xgrid, tt, m),
    col = c("red", "darkgreen", "purple")[m - 2],
    lwd = 2.5,
    lty = 1
  )
}


# Proposed model CDF
lines(
  xgrid,
  pProposed(
    xgrid,
    proposedPars[1],
    proposedPars[2],
    proposedPars[3],
    proposedPars[4]
  ),
  col = "orange",
  lwd = 2.5,
  lty = 1
)


# CDF legend
legend(
  "topleft",
  legend = c(
    "Empirical CDF",
    "Ishita",
    "Polynomial m=3",
    "Polynomial m=4",
    "Polynomial m=5",
    "Proposed MOT2EHL-Gom-Ishita"
  ),
  col = c(
    "black",
    "blue",
    "red",
    "darkgreen",
    "purple",
    "orange"
  ),
  lty = 1,
  lwd = 2.5,
  cex = 0.8,
  bty = "n"
)


# 9. OUTPUT FILE LOCATION


cat("\nOutput directory:\n")
print(getwd())

cat("\nCSV comparison table saved as:\n")
cat("Ishita_model_comparison_geyser.csv\n")








# 1. AML DATA

if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages("survival")
}

library(survival)
data("aml", package = "survival")

time <- as.numeric(aml$time)
status <- as.numeric(aml$status)

keep <- is.finite(time) & time > 0 & status %in% c(0, 1)
time <- time[keep]
status <- status[keep]

# In survival::aml: status = 1 means event; 0 means censored
x <- time
n <- length(x)

cat("Sample size =", n, "\n")
cat("Observed events =", sum(status == 1), "\n")
cat("Censored observations =", sum(status == 0), "\n")
print(summary(x))
print(table(status))


# 2. BASELINE ISHITA PDF AND CDF
dIshita <- function(x, theta) {
  ans <- numeric(length(x))
  ok <- is.finite(x) & x > 0
  
  ans[ok] <- theta^3 / (theta^3 + 2) *
    (theta + x[ok]^2) * exp(-theta * x[ok])
  
  ans
}

pIshita <- function(x, theta) {
  ans <- numeric(length(x))
  ans[x <= 0] <- 0
  ans[is.infinite(x) & x > 0] <- 1
  
  ok <- is.finite(x) & x > 0
  z <- theta * x[ok]
  
  ans[ok] <- theta^3 / (theta^3 + 2) *
    ((1 - exp(-z)) +
       2 / theta^3 *
       (1 - exp(-z) * (1 + z + z^2 / 2)))
  
  pmin(1, pmax(0, ans))
}


# 3. ISHITA-POLYNOMIAL PDF AND CDF

Cm <- function(theta, m) {
  terms <- if (m >= 3) {
    sapply(3:m, function(k) factorial(k) / theta^(k + 1))
  } else {
    numeric(0)
  }
  
  1 / (1 + 2 / theta^3 + sum(terms))
}

dPolyIshita <- function(x, theta, m) {
  ans <- numeric(length(x))
  ok <- is.finite(x) & x > 0
  
  if (any(ok)) {
    xx <- x[ok]
    poly <- theta + xx^2
    
    if (m >= 3) {
      for (k in 3:m) {
        poly <- poly + xx^k
      }
    }
    
    ans[ok] <- Cm(theta, m) * poly * exp(-theta * xx)
  }
  
  ans
}

pPolyIshita <- function(x, theta, m) {
  ans <- numeric(length(x))
  ans[x <= 0] <- 0
  ans[is.infinite(x) & x > 0] <- 1
  
  ok <- is.finite(x) & x > 0
  
  if (any(ok)) {
    z <- theta * x[ok]
    
    val <- 1 - exp(-z)
    val <- val + 2 / theta^3 * pgamma(z, shape = 3)
    
    if (m >= 3) {
      for (k in 3:m) {
        val <- val +
          factorial(k) / theta^(k + 1) *
          pgamma(z, shape = k + 1)
      }
    }
    
    ans[ok] <- Cm(theta, m) * val
  }
  
  pmin(1, pmax(0, ans))
}



# 4. PROPOSED MOT2EHL-GOM-ISHITA MODEL

pProposed <- function(x, delta, alpha, gamma, theta) {
  G <- pIshita(x, theta)
  S <- pmax(1 - G, 1e-14)
  
  v <- (1 - S^(-gamma)) / gamma
  q <- exp(pmax(v, -700))
  A <- (q / (2 - q))^alpha
  
  F <- 1 - delta * A / (1 - (1 - delta) * A)
  
  F[x <= 0] <- 0
  F[is.infinite(x) & x > 0] <- 1
  
  pmin(1, pmax(0, F))
}

dProposed <- function(x, delta, alpha, gamma, theta) {
  G <- pIshita(x, theta)
  g <- dIshita(x, theta)
  S <- pmax(1 - G, 1e-14)
  
  v <- (1 - S^(-gamma)) / gamma
  q <- exp(pmax(v, -700))
  A <- (q / (2 - q))^alpha
  
  den <- 1 - (1 - delta) * A
  
  f <- delta * 2 * alpha * A *
    S^(-gamma - 1) * g /
    ((2 - q) * den^2)
  
  f[!is.finite(f) | x <= 0] <- 0
  pmax(f, 0)
}



# 5. RIGHT-CENSORED MAXIMUM-LIKELIHOOD ESTIMATION
#
# Log-likelihood:
# sum(status*log(f(t)) + (1-status)*log(S(t)))
#
# status = 1: observed event
# status = 0: right-censored

# Ishita model
nll_ishita <- function(par) {
  theta <- exp(par[1])
  
  f <- dIshita(time, theta)
  S <- 1 - pIshita(time, theta)
  
  if (any(!is.finite(f)) || any(f <= 0) ||
      any(!is.finite(S)) || any(S <= 0)) {
    return(1e100)
  }
  
  -sum(status * log(f) + (1 - status) * log(S))
}

fitI <- optim(
  log(1 / mean(time)),
  nll_ishita,
  method = "BFGS",
  hessian = TRUE
)

thetaI <- exp(fitI$par)
llI <- -fitI$value


# Polynomial Ishita models
fitPoly <- function(m) {
  obj <- function(par) {
    theta <- exp(par[1])
    
    f <- dPolyIshita(time, theta, m)
    S <- 1 - pPolyIshita(time, theta, m)
    
    if (any(!is.finite(f)) || any(f <= 0) ||
        any(!is.finite(S)) || any(S <= 0)) {
      return(1e100)
    }
    
    -sum(status * log(f) + (1 - status) * log(S))
  }
  
  fit <- optim(
    log(1 / mean(time)),
    obj,
    method = "BFGS"
  )
  
  list(
    theta = exp(fit$par),
    logLik = -fit$value,
    convergence = fit$convergence
  )
}

polyFits <- lapply(3:5, fitPoly)
names(polyFits) <- paste0("Ishita-Polynomial-m", 3:5)


# Proposed model
nll_proposed <- function(par) {
  pars <- exp(par)
  
  delta <- pars[1]
  alpha <- pars[2]
  gamma <- pars[3]
  theta <- pars[4]
  
  f <- dProposed(time, delta, alpha, gamma, theta)
  S <- 1 - pProposed(time, delta, alpha, gamma, theta)
  
  if (any(!is.finite(f)) || any(f <= 0) ||
      any(!is.finite(S)) || any(S <= 0)) {
    return(1e100)
  }
  
  -sum(status * log(f) + (1 - status) * log(S))
}

fitP <- optim(
  par = log(c(1, 1, 1, 1 / mean(time))),
  fn = nll_proposed,
  method = "Nelder-Mead",
  control = list(maxit = 50000, reltol = 1e-9)
)

proposedPars <- exp(fitP$par)
llP <- -fitP$value


# Display parameter estimates
cat("\nIshita theta estimate:\n")
print(thetaI)

cat("\nPolynomial estimates:\n")
print(lapply(polyFits, function(z) z$theta))

cat("\nProposed estimates (delta, alpha, gamma, theta):\n")
print(proposedPars)

cat("\nProposed model convergence code:", fitP$convergence, "\n")


# 6. MODEL COMPARISON
# LogLik, AIC, BIC, and a descriptive KM-based discrepancy
# Ordinary KS p-values are not calculated for censored data.

modelNames <- c(
  "Ishita",
  names(polyFits),
  "MOT2EHL-Gom-Ishita"
)

logLikValues <- c(
  llI,
  sapply(polyFits, function(z) z$logLik),
  llP
)

numPars <- c(1, 1, 1, 1, 4)

# Fitted CDF functions
cdfFunctions <- c(
  list(function(z) pIshita(z, thetaI)),
  lapply(3:5, function(m) {
    local({
      mm <- m
      tt <- polyFits[[paste0("Ishita-Polynomial-m", mm)]]$theta
      function(z) pPolyIshita(z, tt, mm)
    })
  }),
  list(function(z) {
    pProposed(
      z,
      proposedPars[1],
      proposedPars[2],
      proposedPars[3],
      proposedPars[4]
    )
  })
)

# Kaplan-Meier estimate of the survival function
km <- survfit(Surv(time, status) ~ 1)

# Evaluate fitted CDF and KM CDF at observed times
km_times <- km$time
km_cdf <- 1 - km$surv

# Descriptive discrepancy between KM CDF and model CDF.
# This is NOT the ordinary one-sample KS test statistic.
KM_discrepancy <- sapply(cdfFunctions, function(ff) {
  fitted_cdf <- ff(km_times)
  max(abs(km_cdf - fitted_cdf), na.rm = TRUE)
})

comparison <- data.frame(
  Model = modelNames,
  Parameters = numPars,
  LogLik = logLikValues,
  AIC = -2 * logLikValues + 2 * numPars,
  BIC = -2 * logLikValues + log(n) * numPars,
  KM_discrepancy = KM_discrepancy,
  KS_p_value = NA_real_
)

comparison <- comparison[order(comparison$AIC), ]
rownames(comparison) <- NULL

cat("\nMODEL COMPARISON TABLE\n")
print(comparison, digits = 5)

write.csv(
  comparison,
  "Ishita_model_comparison_aml.csv",
  row.names = FALSE
)


# 7. HISTOGRAM WITH FITTED DENSITIES
# Descriptive only: censored times are not all observed lifetimes.

xgrid <- seq(
  0,
  max(time) * 1.15,
  length.out = 2000
)

dens_I <- dIshita(xgrid, thetaI)

dens_P3 <- dPolyIshita(
  xgrid, polyFits[["Ishita-Polynomial-m3"]]$theta, 3
)

dens_P4 <- dPolyIshita(
  xgrid, polyFits[["Ishita-Polynomial-m4"]]$theta, 4
)

dens_P5 <- dPolyIshita(
  xgrid, polyFits[["Ishita-Polynomial-m5"]]$theta, 5
)

dens_Prop <- dProposed(
  xgrid,
  proposedPars[1],
  proposedPars[2],
  proposedPars[3],
  proposedPars[4]
)

hist_data <- hist(time, breaks = "FD", plot = FALSE)

all_heights <- c(
  hist_data$density, dens_I, dens_P3, dens_P4, dens_P5, dens_Prop
)

all_heights <- all_heights[is.finite(all_heights)]
ylim_max <- max(all_heights, na.rm = TRUE) * 1.15

hist(
  time,
  probability = TRUE,
  breaks = "FD",
  col = "grey90",
  border = "white",
  xlim = c(0, max(time) * 1.15),
  ylim = c(0, ylim_max),
  main = "",
  xlab = "Observed time",
  ylab = "Density",
  cex.main = 1.3,
  cex.lab = 1.2
)

lines(xgrid, dens_I, col = "blue", lwd = 2.5)
lines(xgrid, dens_P3, col = "red", lwd = 2.5)
lines(xgrid, dens_P4, col = "darkgreen", lwd = 2.5)
lines(xgrid, dens_P5, col = "purple", lwd = 2.5)
lines(xgrid, dens_Prop, col = "orange", lwd = 2.5)

legend(
  "topright",
  legend = c(
    "Ishita", "Polynomial m=3", "Polynomial m=4",
    "Polynomial m=5", "Proposed MOT2EHL-Gom-Ishita"
  ),
  col = c("blue", "red", "darkgreen", "purple", "orange"),
  lty = 1,
  lwd = 2.5,
  cex = 0.75,
  bty = "n"
)


# 8. KAPLAN-MEIER SURVIVAL CURVE AND FITTED SURVIVAL CURVES

plot(
  km,
  conf.int = FALSE,
  mark.time = TRUE,
  xlim = c(0, max(time) * 1.15),
  ylim = c(0, 1),
  xlab = "Observed time",
  ylab = "Survival probability",
  main = "",
  lwd = 2,
  col = "black"
)

lines(
  xgrid,
  1 - pIshita(xgrid, thetaI),
  col = "blue",
  lwd = 2.5
)

for (m in 3:5) {
  tt <- polyFits[[paste0("Ishita-Polynomial-m", m)]]$theta
  
  lines(
    xgrid,
    1 - pPolyIshita(xgrid, tt, m),
    col = c("red", "darkgreen", "purple")[m - 2],
    lwd = 2.5
  )
}

lines(
  xgrid,
  1 - pProposed(
    xgrid,
    proposedPars[1],
    proposedPars[2],
    proposedPars[3],
    proposedPars[4]
  ),
  col = "orange",
  lwd = 2.5
)

legend(
  "topright",
  legend = c(
    "Kaplan-Meier",
    "Ishita",
    "Polynomial m=3",
    "Polynomial m=4",
    "Polynomial m=5",
    "Proposed MOT2EHL-Gom-Ishita"
  ),
  col = c(
    "black", "blue", "red", "darkgreen", "purple", "orange"
  ),
  lty = 1,
  lwd = 2.5,
  cex = 0.75,
  bty = "n"
)


# 9. OUTPUT FILE LOCATION

cat("\nOutput directory:\n")
print(getwd())

cat("\nCSV comparison table saved as:\n")
cat("Ishita_model_comparison_aml.csv\n")













# 1. STANFORD2 DATA

if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages("survival")
}

library(survival)
data("stanford2", package = "survival")

time <- as.numeric(stanford2$time)
status <- as.numeric(stanford2$status)

keep <- is.finite(time) & time > 0 & status %in% c(0, 1)
time <- time[keep]
status <- status[keep]

x <- time
n <- length(x)

cat("Sample size =", n, "\n")
cat("Observed events =", sum(status == 1), "\n")
cat("Censored observations =", sum(status == 0), "\n")
print(summary(x))
print(table(status))


# 2. BASELINE ISHITA PDF AND CDF

dIshita <- function(x, theta) {
  ans <- numeric(length(x))
  ok <- is.finite(x) & x > 0
  
  ans[ok] <- theta^3 / (theta^3 + 2) *
    (theta + x[ok]^2) * exp(-theta * x[ok])
  
  ans
}

pIshita <- function(x, theta) {
  ans <- numeric(length(x))
  ans[x <= 0] <- 0
  ans[is.infinite(x) & x > 0] <- 1
  
  ok <- is.finite(x) & x > 0
  z <- theta * x[ok]
  
  ans[ok] <- theta^3 / (theta^3 + 2) *
    ((1 - exp(-z)) +
       2 / theta^3 *
       (1 - exp(-z) * (1 + z + z^2 / 2)))
  
  pmin(1, pmax(0, ans))
}



# 3. ISHITA-POLYNOMIAL PDF AND CDF

Cm <- function(theta, m) {
  terms <- if (m >= 3) {
    sapply(3:m, function(k) factorial(k) / theta^(k + 1))
  } else {
    numeric(0)
  }
  
  1 / (1 + 2 / theta^3 + sum(terms))
}

dPolyIshita <- function(x, theta, m) {
  ans <- numeric(length(x))
  ok <- is.finite(x) & x > 0
  
  if (any(ok)) {
    xx <- x[ok]
    poly <- theta + xx^2
    
    if (m >= 3) {
      for (k in 3:m) {
        poly <- poly + xx^k
      }
    }
    
    ans[ok] <- Cm(theta, m) * poly * exp(-theta * xx)
  }
  
  ans
}

pPolyIshita <- function(x, theta, m) {
  ans <- numeric(length(x))
  ans[x <= 0] <- 0
  ans[is.infinite(x) & x > 0] <- 1
  
  ok <- is.finite(x) & x > 0
  
  if (any(ok)) {
    z <- theta * x[ok]
    
    val <- 1 - exp(-z)
    val <- val + 2 / theta^3 * pgamma(z, shape = 3)
    
    if (m >= 3) {
      for (k in 3:m) {
        val <- val +
          factorial(k) / theta^(k + 1) *
          pgamma(z, shape = k + 1)
      }
    }
    
    ans[ok] <- Cm(theta, m) * val
  }
  
  pmin(1, pmax(0, ans))
}


# 4. PROPOSED MOT2EHL-GOM-ISHITA MODEL

pProposed <- function(x, delta, alpha, gamma, theta) {
  G <- pIshita(x, theta)
  S <- pmax(1 - G, 1e-14)
  
  v <- (1 - S^(-gamma)) / gamma
  q <- exp(pmax(v, -700))
  A <- (q / (2 - q))^alpha
  
  F <- 1 - delta * A / (1 - (1 - delta) * A)
  
  F[x <= 0] <- 0
  F[is.infinite(x) & x > 0] <- 1
  
  pmin(1, pmax(0, F))
}

dProposed <- function(x, delta, alpha, gamma, theta) {
  G <- pIshita(x, theta)
  g <- dIshita(x, theta)
  S <- pmax(1 - G, 1e-14)
  
  v <- (1 - S^(-gamma)) / gamma
  q <- exp(pmax(v, -700))
  A <- (q / (2 - q))^alpha
  
  den <- 1 - (1 - delta) * A
  
  f <- delta * 2 * alpha * A *
    S^(-gamma - 1) * g /
    ((2 - q) * den^2)
  
  f[!is.finite(f) | x <= 0] <- 0
  pmax(f, 0)
}



# 5. RIGHT-CENSORED MAXIMUM-LIKELIHOOD ESTIMATION

# Log-likelihood:
# sum(status*log(f(time)) + (1-status)*log(S(time)))
# status = 1: observed event
# status = 0: right-censored observation

# Fit baseline Ishita model
nll_ishita <- function(par) {
  theta <- exp(par[1])
  
  f <- dIshita(time, theta)
  S <- 1 - pIshita(time, theta)
  
  if (any(!is.finite(f)) || any(f <= 0) ||
      any(!is.finite(S)) || any(S <= 0)) {
    return(1e100)
  }
  
  -sum(status * log(f) + (1 - status) * log(S))
}

fitI <- optim(
  log(1 / mean(time)),
  nll_ishita,
  method = "BFGS",
  hessian = TRUE
)

thetaI <- exp(fitI$par)
llI <- -fitI$value


# Fit polynomial Ishita models for m = 3, 4, 5
fitPoly <- function(m) {
  obj <- function(par) {
    theta <- exp(par[1])
    
    f <- dPolyIshita(time, theta, m)
    S <- 1 - pPolyIshita(time, theta, m)
    
    if (any(!is.finite(f)) || any(f <= 0) ||
        any(!is.finite(S)) || any(S <= 0)) {
      return(1e100)
    }
    
    -sum(status * log(f) + (1 - status) * log(S))
  }
  
  fit <- optim(
    log(1 / mean(time)),
    obj,
    method = "BFGS"
  )
  
  list(
    theta = exp(fit$par),
    logLik = -fit$value,
    convergence = fit$convergence
  )
}

polyFits <- lapply(3:5, fitPoly)
names(polyFits) <- paste0("Ishita-Polynomial-m", 3:5)


# Fit proposed model
nll_proposed <- function(par) {
  pars <- exp(par)
  
  delta <- pars[1]
  alpha <- pars[2]
  gamma <- pars[3]
  theta <- pars[4]
  
  f <- dProposed(time, delta, alpha, gamma, theta)
  S <- 1 - pProposed(time, delta, alpha, gamma, theta)
  
  if (any(!is.finite(f)) || any(f <= 0) ||
      any(!is.finite(S)) || any(S <= 0)) {
    return(1e100)
  }
  
  -sum(status * log(f) + (1 - status) * log(S))
}

fitP <- optim(
  par = log(c(1, 1, 1, 1 / mean(time))),
  fn = nll_proposed,
  method = "Nelder-Mead",
  control = list(maxit = 50000, reltol = 1e-9)
)

proposedPars <- exp(fitP$par)
llP <- -fitP$value


# Display parameter estimates
cat("\nIshita theta estimate:\n")
print(thetaI)

cat("\nPolynomial estimates:\n")
print(lapply(polyFits, function(z) z$theta))

cat("\nProposed estimates (delta, alpha, gamma, theta):\n")
print(proposedPars)

cat("\nProposed model convergence code:", fitP$convergence, "\n")


# 6. MODEL COMPARISON
# LogLik, AIC, BIC, KM discrepancy and KS p-value

modelNames <- c(
  "Ishita",
  names(polyFits),
  "MOT2EHL-Gom-Ishita"
)

logLikValues <- c(
  llI,
  sapply(polyFits, function(z) z$logLik),
  llP
)

numPars <- c(1, 1, 1, 1, 4)


# Fitted CDF functions
cdfFunctions <- c(
  list(function(z) pIshita(z, thetaI)),
  
  lapply(3:5, function(m) {
    local({
      mm <- m
      tt <- polyFits[[paste0("Ishita-Polynomial-m", mm)]]$theta
      function(z) pPolyIshita(z, tt, mm)
    })
  }),
  
  list(function(z) {
    pProposed(
      z,
      proposedPars[1],
      proposedPars[2],
      proposedPars[3],
      proposedPars[4]
    )
  })
)


# Kaplan-Meier estimate
km <- survfit(Surv(time, status) ~ 1)

km_times <- km$time
km_cdf <- 1 - km$surv


# Descriptive discrepancy between KM CDF and fitted CDF
# This is not the ordinary one-sample KS statistic.
KM_discrepancy <- sapply(cdfFunctions, function(ff) {
  fitted_cdf <- ff(km_times)
  max(abs(km_cdf - fitted_cdf), na.rm = TRUE)
})


# Ordinary KS p-values are not valid for censored observations.
# Use NA rather than reporting misleading p-values.
KS_p_value <- rep(NA_real_, length(modelNames))

comparison <- data.frame(
  Model = modelNames,
  Parameters = numPars,
  LogLik = logLikValues,
  AIC = -2 * logLikValues + 2 * numPars,
  BIC = -2 * logLikValues + log(n) * numPars,
  KM_discrepancy = KM_discrepancy,
  KS_p_value = KS_p_value
)

comparison <- comparison[order(comparison$AIC), ]
rownames(comparison) <- NULL

cat("\nMODEL COMPARISON TABLE\n")
print(comparison, digits = 5)

write.csv(
  comparison,
  "Ishita_model_comparison_stanford2.csv",
  row.names = FALSE
)

# 7. HISTOGRAM WITH FITTED DENSITIES
# Descriptive plot of observed follow-up times

xgrid <- seq(
  0,
  max(time) * 1.15,
  length.out = 2000
)

dens_I <- dIshita(xgrid, thetaI)

dens_P3 <- dPolyIshita(
  xgrid,
  polyFits[["Ishita-Polynomial-m3"]]$theta,
  3
)

dens_P4 <- dPolyIshita(
  xgrid,
  polyFits[["Ishita-Polynomial-m4"]]$theta,
  4
)

dens_P5 <- dPolyIshita(
  xgrid,
  polyFits[["Ishita-Polynomial-m5"]]$theta,
  5
)

dens_Prop <- dProposed(
  xgrid,
  proposedPars[1],
  proposedPars[2],
  proposedPars[3],
  proposedPars[4]
)

hist_data <- hist(
  time,
  breaks = "FD",
  plot = FALSE
)

all_heights <- c(
  hist_data$density,
  dens_I,
  dens_P3,
  dens_P4,
  dens_P5,
  dens_Prop
)

all_heights <- all_heights[is.finite(all_heights)]
ylim_max <- max(all_heights, na.rm = TRUE) * 1.15

hist(
  time,
  probability = TRUE,
  breaks = "FD",
  col = "grey90",
  border = "white",
  xlim = c(0, max(time) * 1.15),
  ylim = c(0, ylim_max),
  main = "",
  xlab = "Survival or censoring time",
  ylab = "Density",
  cex.main = 1.3,
  cex.lab = 1.2
)

lines(xgrid, dens_I, col = "blue", lwd = 2.5)
lines(xgrid, dens_P3, col = "red", lwd = 2.5)
lines(xgrid, dens_P4, col = "darkgreen", lwd = 2.5)
lines(xgrid, dens_P5, col = "purple", lwd = 2.5)
lines(xgrid, dens_Prop, col = "orange", lwd = 2.5)

legend(
  "topright",
  legend = c(
    "Ishita",
    "Polynomial m=3",
    "Polynomial m=4",
    "Polynomial m=5",
    "Proposed MOT2EHL-Gom-Ishita"
  ),
  col = c("blue", "red", "darkgreen", "purple", "orange"),
  lty = 1,
  lwd = 2.5,
  cex = 0.75,
  bty = "n"
)


# 8. KAPLAN-MEIER AND FITTED SURVIVAL CURVES

plot(
  km,
  conf.int = FALSE,
  mark.time = TRUE,
  xlim = c(0, max(time) * 1.15),
  ylim = c(0, 1),
  xlab = "Survival or censoring time",
  ylab = "Survival probability",
  main = "",
  lwd = 2,
  col = "black"
)

lines(
  xgrid,
  1 - pIshita(xgrid, thetaI),
  col = "blue",
  lwd = 2.5
)

for (m in 3:5) {
  tt <- polyFits[[paste0("Ishita-Polynomial-m", m)]]$theta
  
  lines(
    xgrid,
    1 - pPolyIshita(xgrid, tt, m),
    col = c("red", "darkgreen", "purple")[m - 2],
    lwd = 2.5
  )
}

lines(
  xgrid,
  1 - pProposed(
    xgrid,
    proposedPars[1],
    proposedPars[2],
    proposedPars[3],
    proposedPars[4]
  ),
  col = "orange",
  lwd = 2.5
)

legend(
  "topright",
  legend = c(
    "Kaplan-Meier",
    "Ishita",
    "Polynomial m=3",
    "Polynomial m=4",
    "Polynomial m=5",
    "Proposed MOT2EHL-Gom-Ishita"
  ),
  col = c(
    "black", "blue", "red", "darkgreen", "purple", "orange"
  ),
  lty = 1,
  lwd = 2.5,
  cex = 0.75,
  bty = "n"
)


# 9. OUTPUT FILE LOCATION

cat("\nOutput directory:\n")
print(getwd())

cat("\nCSV comparison table saved as:\n")
cat("Ishita_model_comparison_stanford2.csv\n")














# ISHITA AND ISHITA-POLYNOMIAL DISTRIBUTIONS
# PROPOSED MOT2EHL-GOM-ISHITA DISTRIBUTION
# APPLICATION TO THE BUILT-IN LUNG LIFETIME DATASET
#
# Dataset: survival::lung
# Original status: 1 = censored, 2 = dead
# Converted status: 0 = censored, 1 = event
#
# Methods:
#   1. Ishita distribution
#   2. Ishita-polynomial distribution, m = 3
#   3. Ishita-polynomial distribution, m = 4
#   4. Ishita-polynomial distribution, m = 5
#   5. Proposed MOT2EHL-Gom-Ishita distribution
#
# Estimation: Maximum likelihood with right censoring


# 1. INSTALL AND LOAD PACKAGES

if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages("survival")
}

library(survival)


# 2. LOAD AND PREPARE THE LUNG DATASET

data("lung", package = "survival")

time <- as.numeric(lung$time)

# Convert original coding:
# 1 = censored, 2 = dead
# New coding:
# 0 = censored, 1 = event

event <- as.integer(lung$status == 2)

# Remove missing or invalid observations

valid <- is.finite(time) & time > 0 & !is.na(event)

time <- time[valid]
event <- event[valid]

n <- length(time)

cat("\n========================================\n")
cat("LUNG DATASET SUMMARY\n")
cat("========================================\n")
cat("Number of observations:", n, "\n")
cat("Number of observed deaths:", sum(event == 1), "\n")
cat("Number of censored observations:", sum(event == 0), "\n")
cat("Minimum follow-up time:", min(time), "days\n")
cat("Maximum follow-up time:", max(time), "days\n")
cat("Mean follow-up time:", mean(time), "days\n")
cat("========================================\n\n")



# 3. ISHITA DISTRIBUTION

# PDF:
# f(x) = theta^3/(theta^3+2) *
#        (theta+x^2)*exp(-theta*x), x > 0

dIshita <- function(x, theta) {
  
  ans <- numeric(length(x))
  
  ok <- is.finite(x) & x > 0
  
  ans[ok] <-
    theta^3 / (theta^3 + 2) *
    (theta + x[ok]^2) *
    exp(-theta * x[ok])
  
  ans
}


# CDF of Ishita distribution

pIshita <- function(x, theta) {
  
  ans <- numeric(length(x))
  
  ans[x <= 0] <- 0
  ans[is.infinite(x) & x > 0] <- 1
  
  ok <- is.finite(x) & x > 0
  
  if (any(ok)) {
    
    z <- theta * x[ok]
    
    ans[ok] <-
      theta^3 / (theta^3 + 2) *
      (
        1 - exp(-z) +
          2 / theta^3 *
          (1 - exp(-z) *
             (1 + z + z^2 / 2))
      )
  }
  
  pmin(1, pmax(0, ans))
}



# 4. ISHITA-POLYNOMIAL DISTRIBUTION

# PDF:
# f_m(x) = C_m(theta) *
#          (theta+x^2+sum(x^k, k=3,...,m)) *
#          exp(-theta*x)
#
# Normalizing constant:
# C_m(theta) = 1 /
# [1 + 2/theta^3 + sum(k!/theta^(k+1))]


Cm <- function(theta, m) {
  
  terms <- if (m >= 3) {
    sapply(3:m, function(k) {
      factorial(k) / theta^(k + 1)
    })
  } else {
    numeric(0)
  }
  
  1 / (1 + 2 / theta^3 + sum(terms))
}


# Ishita-polynomial PDF

dPolyIshita <- function(x, theta, m) {
  
  ans <- numeric(length(x))
  
  ok <- is.finite(x) & x > 0
  
  if (any(ok)) {
    
    xx <- x[ok]
    
    poly_part <- theta + xx^2
    
    if (m >= 3) {
      for (k in 3:m) {
        poly_part <- poly_part + xx^k
      }
    }
    
    ans[ok] <-
      Cm(theta, m) *
      poly_part *
      exp(-theta * xx)
  }
  
  ans
}


# Ishita-polynomial CDF

pPolyIshita <- function(x, theta, m) {
  
  ans <- numeric(length(x))
  
  ans[x <= 0] <- 0
  ans[is.infinite(x) & x > 0] <- 1
  
  ok <- is.finite(x) & x > 0
  
  if (any(ok)) {
    
    xx <- x[ok]
    z <- theta * xx
    
    cdf_part <-
      1 - exp(-z) +
      2 / theta^3 * pgamma(z, shape = 3)
    
    if (m >= 3) {
      
      for (k in 3:m) {
        
        cdf_part <- cdf_part +
          factorial(k) / theta^(k + 1) *
          pgamma(z, shape = k + 1)
      }
    }
    
    ans[ok] <- Cm(theta, m) * cdf_part
  }
  
  pmin(1, pmax(0, ans))
}



# 5. PROPOSED MOT2EHL-GOM-ISHITA DISTRIBUTION

# Parameters:
# delta > 0, alpha > 0, gamma > 0, theta > 0
#
# The following uses the proposed CDF/PDF formulation from
# the previous model-comparison code.


pProposed <- function(x, delta, alpha, gamma, theta) {
  
  G <- pIshita(x, theta)
  
  S <- pmax(1 - G, 1e-14)
  
  v <- (1 - S^(-gamma)) / gamma
  
  q <- exp(pmax(v, -700))
  
  A <- (q / (2 - q))^alpha
  
  den <- 1 - (1 - delta) * A
  
  F <- 1 - delta * A / den
  
  F[x <= 0] <- 0
  
  F[is.infinite(x) & x > 0] <- 1
  
  pmin(1, pmax(0, F))
}


# Proposed PDF

dProposed <- function(x, delta, alpha, gamma, theta) {
  
  G <- pIshita(x, theta)
  
  g <- dIshita(x, theta)
  
  S <- pmax(1 - G, 1e-14)
  
  v <- (1 - S^(-gamma)) / gamma
  
  q <- exp(pmax(v, -700))
  
  A <- (q / (2 - q))^alpha
  
  den <- 1 - (1 - delta) * A
  
  f <- delta * 2 * alpha * A *
    S^(-gamma - 1) * g /
    ((2 - q) * den^2)
  
  f[!is.finite(f) | x <= 0] <- 0
  
  pmax(f, 0)
}


# 6. RIGHT-CENSORED LOG-LIKELIHOOD

# Log-likelihood:
#
# L = sum[event*log(f(time)) +
#         (1-event)*log(S(time))]
#
# where S(time) = 1 - F(time)


censored_loglik <- function(par, model, time, event) {
  
  # Parameter order for each model:
  # Ishita: theta
  # Polynomial: theta
  # Proposed: delta, alpha, gamma, theta
  
  if (any(!is.finite(par)) || any(par <= 0)) {
    return(-1e100)
  }
  
  if (model == "Ishita") {
    
    theta <- par[1]
    
    f <- dIshita(time, theta)
    F <- pIshita(time, theta)
    
  } else if (model == "Poly3") {
    
    theta <- par[1]
    
    f <- dPolyIshita(time, theta, 3)
    F <- pPolyIshita(time, theta, 3)
    
  } else if (model == "Poly4") {
    
    theta <- par[1]
    
    f <- dPolyIshita(time, theta, 4)
    F <- pPolyIshita(time, theta, 4)
    
  } else if (model == "Poly5") {
    
    theta <- par[1]
    
    f <- dPolyIshita(time, theta, 5)
    F <- pPolyIshita(time, theta, 5)
    
  } else if (model == "Proposed") {
    
    delta <- par[1]
    alpha <- par[2]
    gamma <- par[3]
    theta <- par[4]
    
    f <- dProposed(
      time, delta, alpha, gamma, theta
    )
    
    F <- pProposed(
      time, delta, alpha, gamma, theta
    )
    
  } else {
    stop("Unknown model")
  }
  
  S <- 1 - F
  
  # Protect against log(0) and numerical underflow
  
  f <- pmax(f, 1e-300)
  S <- pmax(S, 1e-12)
  
  loglik <- sum(
    event * log(f) +
      (1 - event) * log(S)
  )
  
  if (!is.finite(loglik)) {
    return(-1e100)
  }
  
  loglik
}


# 7. MAXIMUM LIKELIHOOD ESTIMATION

fit_model <- function(model, time, event) {
  
  # Starting values and parameter bounds
  
  theta_start <- 1 / mean(time)
  
  if (model == "Ishita") {
    
    start <- c(theta_start)
    lower <- c(1e-6)
    upper <- c(10)
    
  } else if (model %in% c("Poly3", "Poly4", "Poly5")) {
    
    start <- c(theta_start)
    lower <- c(1e-6)
    upper <- c(10)
    
  } else if (model == "Proposed") {
    
    start <- c(
      delta = 1,
      alpha = 1,
      gamma = 1,
      theta = theta_start
    )
    
    lower <- c(1e-4, 1e-4, 1e-4, 1e-6)
    upper <- c(100, 100, 100, 10)
  }
  
  # Optimize on log-parameter scale to enforce positivity
  
  objective <- function(logpar) {
    
    par <- exp(logpar)
    
    ll <- censored_loglik(
      par = par,
      model = model,
      time = time,
      event = event
    )
    
    -ll
  }
  
  opt <- optim(
    par = log(start),
    fn = objective,
    method = "L-BFGS-B",
    lower = log(lower),
    upper = log(upper),
    control = list(
      maxit = 5000,
      factr = 1e7
    )
  )
  
  estimates <- exp(opt$par)
  
  names(estimates) <- if (model == "Proposed") {
    c("delta", "alpha", "gamma", "theta")
  } else {
    "theta"
  }
  
  ll <- censored_loglik(
    par = estimates,
    model = model,
    time = time,
    event = event
  )
  
  k <- length(estimates)
  
  AIC_value <- -2 * ll + 2 * k
  
  BIC_value <- -2 * ll + log(length(time)) * k
  
  list(
    model = model,
    estimates = estimates,
    logLik = ll,
    AIC = AIC_value,
    BIC = BIC_value,
    convergence = opt$convergence,
    message = opt$message
  )
}


# 8. FIT ALL FIVE MODELS

models <- c(
  "Ishita",
  "Poly3",
  "Poly4",
  "Poly5",
  "Proposed"
)

fits <- lapply(
  models,
  function(mod) {
    fit_model(
      model = mod,
      time = time,
      event = event
    )
  }
)

names(fits) <- models


# 9. DISPLAY PARAMETER ESTIMATES

cat("\n========================================\n")
cat("MAXIMUM LIKELIHOOD ESTIMATES\n")
cat("========================================\n")

for (mod in models) {
  
  cat("\nModel:", mod, "\n")
  
  print(fits[[mod]]$estimates)
  
  cat("Log-likelihood:", fits[[mod]]$logLik, "\n")
  
  cat("Convergence code:", fits[[mod]]$convergence, "\n")
}



# 10. MODEL COMPARISON TABLE

parameter_count <- c(
  Ishita = 1,
  Poly3 = 1,
  Poly4 = 1,
  Poly5 = 1,
  Proposed = 4
)

comparison <- data.frame(
  Model = models,
  Parameters = as.integer(parameter_count[models]),
  LogLik = sapply(fits, function(z) z$logLik),
  AIC = sapply(fits, function(z) z$AIC),
  BIC = sapply(fits, function(z) z$BIC),
  Convergence = sapply(fits, function(z) z$convergence),
  stringsAsFactors = FALSE
)

comparison <- comparison[
  order(comparison$AIC),
]

rownames(comparison) <- NULL

cat("\n========================================\n")
cat("MODEL COMPARISON TABLE\n")
cat("========================================\n")

print(comparison, digits = 6)

# Save results to CSV

write.csv(
  comparison,
  "Ishita_model_comparison_lung.csv",
  row.names = FALSE
)

cat("\nComparison table saved as:\n")
cat("Ishita_model_comparison_lung.csv\n")


# 11. CREATE CDF AND PDF FUNCTIONS FOR EACH FITTED MODEL

cdfFunctions <- list(
  
  Ishita = function(x) {
    pIshita(x, fits$Ishita$estimates[1])
  },
  
  Poly3 = function(x) {
    pPolyIshita(x, fits$Poly3$estimates[1], 3)
  },
  
  Poly4 = function(x) {
    pPolyIshita(x, fits$Poly4$estimates[1], 4)
  },
  
  Poly5 = function(x) {
    pPolyIshita(x, fits$Poly5$estimates[1], 5)
  },
  
  Proposed = function(x) {
    pProposed(
      x,
      fits$Proposed$estimates[1],
      fits$Proposed$estimates[2],
      fits$Proposed$estimates[3],
      fits$Proposed$estimates[4]
    )
  }
)


pdfFunctions <- list(
  
  Ishita = function(x) {
    dIshita(x, fits$Ishita$estimates[1])
  },
  
  Poly3 = function(x) {
    dPolyIshita(x, fits$Poly3$estimates[1], 3)
  },
  
  Poly4 = function(x) {
    dPolyIshita(x, fits$Poly4$estimates[1], 4)
  },
  
  Poly5 = function(x) {
    dPolyIshita(x, fits$Poly5$estimates[1], 5)
  },
  
  Proposed = function(x) {
    dProposed(
      x,
      fits$Proposed$estimates[1],
      fits$Proposed$estimates[2],
      fits$Proposed$estimates[3],
      fits$Proposed$estimates[4]
    )
  }
)



# 12. KAPLAN-MEIER ESTIMATE

km_fit <- survfit(
  Surv(time, event) ~ 1
)

km_times <- km_fit$time

km_cdf <- 1 - km_fit$surv


# Compare model CDFs with the Kaplan-Meier CDF

KM_discrepancy <- sapply(
  cdfFunctions,
  function(ff) {
    max(
      abs(km_cdf - ff(km_times)),
      na.rm = TRUE
    )
  }
)

comparison$KM_discrepancy <- KM_discrepancy[
  comparison$Model
]

# Save updated table

write.csv(
  comparison,
  "Ishita_model_comparison_lung.csv",
  row.names = FALSE
)

cat("\n========================================\n")
cat("MODEL COMPARISON WITH KM DISCREPANCY\n")
cat("========================================\n")

print(comparison, digits = 6)


# 13. HISTOGRAM WITH FITTED DENSITIES

# CAUTION:
# Observed follow-up times include censored observations.
# This histogram is descriptive and is not a histogram of
# fully observed death times.

x_grid <- seq(
  0.01,
  max(time),
  length.out = 1000
)

hist(
  time,
  breaks = 25,
  probability = TRUE,
  col = "grey85",
  border = "white",
  main = "",
  xlab = "Observed follow-up time (days)",
  ylab = "Density"
)

line_colors <- c(
  "blue",
  "red",
  "darkgreen",
  "purple",
  "orange"
)

for (i in seq_along(models)) {
  
  mod <- models[i]
  
  lines(
    x_grid,
    pdfFunctions[[mod]](x_grid),
    col = line_colors[i],
    lwd = 2
  )
}

legend(
  "topright",
  legend = models,
  col = line_colors,
  lwd = 2,
  cex = 0.8,
  bty = "n"
)


# 14. KAPLAN-MEIER SURVIVAL CURVE WITH FITTED MODELS

plot(
  km_fit,
  conf.int = FALSE,
  mark.time = TRUE,
  lwd = 2,
  col = "black",
  xlab = "Survival time (days)",
  ylab = "Survival probability",
  main = ""
)

# Draw fitted survival functions: S(t) = 1 - F(t)

for (i in seq_along(models)) {
  
  mod <- models[i]
  
  fitted_survival <- 1 - cdfFunctions[[mod]](x_grid)
  
  lines(
    x_grid,
    pmax(0, pmin(1, fitted_survival)),
    col = line_colors[i],
    lwd = 2
  )
}

legend(
  "topright",
  legend = c("Kaplan-Meier", models),
  col = c("black", line_colors),
  lwd = 2,
  cex = 0.75,
  bty = "n"
)


# 15. EMPIRICAL CDF AND FITTED CDFS


plot(
  km_times,
  km_cdf,
  type = "s",
  lwd = 2,
  col = "black",
  xlab = "Survival time (days)",
  ylab = "Cumulative probability",
  main = "",
  ylim = c(0, 1)
)

for (i in seq_along(models)) {
  
  mod <- models[i]
  
  lines(
    x_grid,
    cdfFunctions[[mod]](x_grid),
    col = line_colors[i],
    lwd = 2
  )
}

legend(
  "bottomright",
  legend = c("Kaplan-Meier CDF", models),
  col = c("black", line_colors),
  lwd = 2,
  cex = 0.75,
  bty = "n"
)


# 16. SAVE MODEL ESTIMATES TO A CSV FILE


estimate_rows <- lapply(
  models,
  function(mod) {
    
    est <- fits[[mod]]$estimates
    
    data.frame(
      Model = mod,
      Parameter = names(est),
      Estimate = as.numeric(est),
      row.names = NULL
    )
  }
)

estimates_table <- do.call(
  rbind,
  estimate_rows
)

write.csv(
  estimates_table,
  "Ishita_parameter_estimates_lung.csv",
  row.names = FALSE
)

cat("\n========================================\n")
cat("ANALYSIS COMPLETED\n")
cat("========================================\n")
cat("Files saved:\n")
cat("1. Ishita_model_comparison_lung.csv\n")
cat("2. Ishita_parameter_estimates_lung.csv\n")
cat("\nPlots are displayed in the R graphics window.\n")

