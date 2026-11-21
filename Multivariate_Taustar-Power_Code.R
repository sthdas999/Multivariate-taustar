
# ============================================================
# POWER ANALYSIS OF THE MULTIVARIATE TAU-STAR TEST
#
# Two asymptotic regimes:
#
# CASE I:
# First-order nondegenerate U-statistic
# sqrt(n) U_n -> N(mu, sigma1^2)
#
# CASE II:
# Second-order degenerate U-statistic
# n U_n -> 6 sum lambda_r (Z_r^2 - 1)
#
# Under contiguous alternatives in Case II:
# Z_r ~ N(delta * c_r, 1)
#
# IMPORTANT:
# The eigenvalues and score coefficients below are illustrative.
# They must eventually be estimated/derived from the actual
# second-order Hoeffding projection of the proposed kernel.
# ============================================================

set.seed(123)

# ------------------------------------------------------------
# GENERAL PARAMETERS
# ------------------------------------------------------------

alpha <- 0.05
B <- 100000

# Choose the asymptotic regime:
# "first_order" or "second_order"

analysis_case <- "second_order"

# ============================================================
# CASE I: FIRST-ORDER NONDEGENERATE REGIME
# ============================================================

# Under H0:
# sqrt(n) U_n -> N(0, sigma1^2)
#
# Under contiguous alternatives:
# sqrt(n) U_n -> N(mu, sigma1^2)

sigma1 <- 1

mu_values <- seq(0, 4, by = 0.5)

# Two-sided critical value
z_critical <- qnorm(1 - alpha / 2)

# ------------------------------------------------------------
# 1A. Theoretical power
# ------------------------------------------------------------

theoretical_power_first <- function(mu, sigma, alpha) {

  z <- qnorm(1 - alpha / 2)

  1 - pnorm(z, mean = mu / sigma) +
    pnorm(-z, mean = mu / sigma)

}

power_first <- sapply(
  mu_values,
  theoretical_power_first,
  sigma = sigma1,
  alpha = alpha
)

# ------------------------------------------------------------
# 1B. Monte Carlo power
# ------------------------------------------------------------

simulated_power_first <- sapply(mu_values, function(mu) {

  test_statistics <- rnorm(
    n = B,
    mean = mu,
    sd = sigma1
  )

  mean(abs(test_statistics) > z_critical * sigma1)

})

# ------------------------------------------------------------
# 1C. First-order power table
# ------------------------------------------------------------

results_first <- data.frame(
  Mean_Shift = mu_values,
  Theoretical_Power = power_first,
  Simulated_Power = simulated_power_first
)

results_first$Theoretical_Power <-
  round(results_first$Theoretical_Power, 4)

results_first$Simulated_Power <-
  round(results_first$Simulated_Power, 4)

cat("\nFIRST-ORDER NONDEGENERATE CASE\n")
print(results_first, row.names = FALSE)


# ============================================================
# CASE II: SECOND-ORDER DEGENERATE REGIME
# ============================================================

# Under H0:
#
# n U_n -> 6 sum lambda_r (Z_r^2 - 1)
#
# Under contiguous alternatives:
#
# n U_n -> 6 sum lambda_r
#              [(Z_r + delta*c_r)^2 - 1]
#
# where Z_r are independent N(0,1).
#
# The eigenvalues and score coefficients below are
# illustrative and must be replaced by theoretically
# justified values for the proposed kernel.

# ------------------------------------------------------------
# 2A. Eigenvalues of the second-order projection
# ------------------------------------------------------------

lambda_values <- c(0.60, 0.25, 0.10, 0.05)

# Corresponding score coefficients:
# c_r = E[s(Z) phi_r(Z)]

score_coefficients <- c(
  0.70, 0.45, 0.25, 0.15
)

# Check square integrability of the score coefficients
if (sum(score_coefficients^2) > 1) {
  stop("The squared score coefficients must not exceed 1.")
}

# Check eigenvalue specifications
if (any(lambda_values < 0)) {
  stop("This implementation assumes nonnegative eigenvalues.")
}

if (length(lambda_values) != length(score_coefficients)) {
  stop("Eigenvalues and score coefficients must have equal length.")
}

R <- length(lambda_values)

# Local alternative parameter
delta_values <- seq(0, 4, by = 0.5)

# ------------------------------------------------------------
# 2B. Generate the limiting second-order statistic
# ------------------------------------------------------------

generate_second_order <- function(B,
                                  delta,
                                  lambda,
                                  score_coefficients) {

  R <- length(lambda)

  Z <- matrix(
    rnorm(B * R),
    nrow = B,
    ncol = R
  )

  # Mean shift under contiguous alternatives
  Z_alt <- sweep(
    Z,
    MARGIN = 2,
    STATS = delta * score_coefficients,
    FUN = "+"
  )

  # Limiting statistic:
  # 6 sum lambda_r [(Z_r + delta*c_r)^2 - 1]

  statistic <- 6 * rowSums(
    sweep(
      Z_alt^2 - 1,
      MARGIN = 2,
      STATS = lambda,
      FUN = "*"
    )
  )

  statistic
}

# ------------------------------------------------------------
# 2C. Null distribution and critical value
# ------------------------------------------------------------

null_statistics <- generate_second_order(
  B = B,
  delta = 0,
  lambda = lambda_values,
  score_coefficients = score_coefficients
)

# Two-sided critical value based on the null distribution
critical_value <- quantile(
  abs(null_statistics),
  probs = 1 - alpha,
  names = FALSE
)

cat("\nSECOND-ORDER DEGENERATE CASE\n")
cat("Two-sided critical value =", critical_value, "\n")

# ------------------------------------------------------------
# 2D. Monte Carlo power under contiguous alternatives
# ------------------------------------------------------------

simulated_power_second <- sapply(delta_values, function(delta) {

  alternative_statistics <- generate_second_order(
    B = B,
    delta = delta,
    lambda = lambda_values,
    score_coefficients = score_coefficients
  )

  mean(abs(alternative_statistics) > critical_value)

})

# ------------------------------------------------------------
# 2E. Second-order power table
# ------------------------------------------------------------

results_second <- data.frame(
  Local_Parameter = delta_values,
  Simulated_Power = simulated_power_second
)

results_second$Simulated_Power <-
  round(results_second$Simulated_Power, 4)

print(results_second, row.names = FALSE)


# ============================================================
# 3. FIRST-ORDER THEORETICAL POWER CURVE
# ============================================================

mu_grid <- seq(0, 4, length.out = 500)

power_grid_first <- sapply(
  mu_grid,
  theoretical_power_first,
  sigma = sigma1,
  alpha = alpha
)

jpeg(
  filename = "power-curve-taustar-first-order.jpg",
  width = 1800,
  height = 1200,
  res = 200,
  quality = 100
)

plot(
  mu_grid,
  power_grid_first,
  type = "l",
  lwd = 2,
  xlim = c(0, 4),
  ylim = c(0, 1),
  xlab = expression(mu),
  ylab = "Power",
  main = expression(
    "First-order asymptotic power of " ~ hat(tau)^"*"
  )
)

abline(h = alpha, lty = 2)

points(
  mu_values,
  power_first,
  pch = 16
)

legend(
  "bottomright",
  legend = c(
    "Theoretical power",
    "Significance level"
  ),
  lty = c(1, 2),
  pch = c(16, NA),
  lwd = c(2, 1),
  bty = "n"
)

dev.off()


# ============================================================
# 4. SECOND-ORDER POWER CURVE
# ============================================================

delta_grid <- seq(0, 4, length.out = 500)

# Use the same null reference distribution and critical value
power_grid_second <- sapply(delta_grid, function(delta) {

  alternative_statistics <- generate_second_order(
    B = B,
    delta = delta,
    lambda = lambda_values,
    score_coefficients = score_coefficients
  )

  mean(abs(alternative_statistics) > critical_value)

})

jpeg(
  filename = "power-curve-taustar-second-order.jpg",
  width = 1800,
  height = 1200,
  res = 200,
  quality = 100
)

plot(
  delta_grid,
  power_grid_second,
  type = "l",
  lwd = 2,
  xlim = c(0, 4),
  ylim = c(0, 1),
  xlab = expression(delta),
  ylab = "Power",
  main = expression(
    "Second-order asymptotic power of " ~ hat(tau)^"*"
  )
)

abline(h = alpha, lty = 2)

points(
  delta_values,
  simulated_power_second,
  pch = 16
)

legend(
  "bottomright",
  legend = c(
    "Asymptotic Monte Carlo power",
    "Significance level"
  ),
  lty = c(1, 2),
  pch = c(16, NA),
  lwd = c(2, 1),
  bty = "n"
)

dev.off()


# ============================================================
# 5. SELECTED ANALYSIS OUTPUT
# ============================================================

if (analysis_case == "first_order") {

  cat("\nSelected regime: FIRST-ORDER NONDEGENERATE\n")
  print(results_first, row.names = FALSE)

} else if (analysis_case == "second_order") {

  cat("\nSelected regime: SECOND-ORDER DEGENERATE\n")
  print(results_second, row.names = FALSE)

} else {

  stop("Choose either 'first_order' or 'second_order'.")

}

cat("\nPower analysis completed.\n")
