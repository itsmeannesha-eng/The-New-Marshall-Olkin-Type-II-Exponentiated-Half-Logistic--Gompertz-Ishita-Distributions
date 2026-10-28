# Baseline Ishita PDF and CDF
g_ishita <- function(x, theta) {
  ifelse(x <= 0, 0, (theta^3 / (theta^3 + 2)) * (theta + x^2) * exp(-theta * x))
}

G_ishita <- function(x, theta) {
  ifelse(x <= 0, 0, 1 - (1 + (theta * x * (theta * x + 2)) / (theta^3 + 2)) * exp(-theta * x))
}

# Vectorized HRF Function for MO-T2EHL-Gom-Ishita
hrf_MOT2EHL_Gom_Ishita <- function(x, delta, alpha, gamma, theta) {
  res <- numeric(length(x))
  pos <- x > 0
  if (!any(pos)) return(res)
  
  x_pos <- x[pos]
  d_bar <- 1 - delta
  
  G_val <- G_ishita(x_pos, theta)
  g_val <- g_ishita(x_pos, theta)
  
  S_bar <- 1 - G_val
  
  # Inner Transformations
  A_x <- 1 - (S_bar)^(-gamma)
  B_x <- exp((1 / gamma) * A_x)
  C_x <- B_x / (2 - B_x)
  
  # HRF Assembly (using the corrected -1 exponent term)
  term1 <- 2 * alpha
  term2 <- exp((alpha / gamma) * A_x)
  term3 <- S_bar^(-gamma - 1)
  term4 <- (2 - B_x)^(-(alpha + 1))
  term5 <- g_val / (C_x^alpha)
  term6 <- (1 - d_bar * (C_x)^alpha)^(-1)
  
  res[pos] <- term1 * term2 * term3 * term4 * term5 * term6
  return(res)
}

# Parameter Vectors (Same 5 Parameter Sets)
delta <- c(1.2, 0.5, 2.5, 0.8, 1.0)
alpha <- c(1.5, 0.4, 2.5, 0.7, 1.0)
gamma <- c(0.8, 0.5, 0.8, 1.5, 1.0)
theta <- c(1.0, 2.0, 1.0, 0.6, 1.2)

# Evaluation sequence
x_seq <- seq(0.01, 3.5, length.out = 600)
colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")

# Graph Setup
plot(NULL, xlim = c(0, 3.5), ylim = c(0, 3), 
     xlab = "x", ylab = "h(x)", 
     main = "",
     cex.main = 1.2, font.main = 2)
grid(col = "gray90", lty = "solid")

# Curve generation
for (i in 1:5) {
  y_vals <- hrf_MOT2EHL_Gom_Ishita(x_seq, delta[i], alpha[i], gamma[i], theta[i])
  lines(x_seq, y_vals, col = colors[i], lwd = 2.5)
}

# Legend
legend("topright", 
       legend = c(
         expression(paste("", delta==1.2, ", ", alpha==1.5, ", ", gamma==0.8, ", ", theta==1.0)),
         expression(paste("", delta==0.5, ", ", alpha==0.4, ", ", gamma==0.5, ", ", theta==2.0)),
         expression(paste("", delta==2.5, ", ", alpha==2.5, ", ", gamma==0.8, ", ", theta==1.0)),
         expression(paste("", delta==0.8, ", ", alpha==0.7, ", ", gamma==1.5, ", ", theta==0.6)),
         expression(paste("", delta==1.0, ", ", alpha==1.0, ", ", gamma==1.0, ", ", theta==1.2))
       ),
       col = colors, lwd = 2.5, bty = "n", cex = 0.85)

