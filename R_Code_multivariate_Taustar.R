# ============================================================
# Power analysis of the multivariate tau-star independence test
# Theoretical and Monte Carlo power under contiguous alternatives
# ============================================================

set.seed(123)

# Parameters
alpha <- 0.05
sigma <- 1
B <- 100000
mu_values <- seq(0, 4, by = 0.5)

# One-sided critical value
z_critical <- qnorm(1 - alpha)

# ------------------------------------------------------------
# 1. Theoretical power
# ------------------------------------------------------------

theoretical_power <- pnorm(
  -z_critical + mu_values / sigma
)

# ------------------------------------------------------------
# 2. Monte Carlo simulation of power
# ------------------------------------------------------------

simulated_power <- sapply(mu_values, function(mu) {
  test_statistics <- rnorm(
    n = B,
    mean = mu / sigma,
    sd = 1
  )
  
  mean(test_statistics > z_critical)
})

# Numerical power table
power_results <- data.frame(
  Mu = mu_values,
  Theoretical_Power = theoretical_power,
  Simulated_Power = simulated_power
)

power_results$Theoretical_Power <-
  round(power_results$Theoretical_Power, 3)

power_results$Simulated_Power <-
  round(power_results$Simulated_Power, 3)

print(power_results, row.names = FALSE)

# ------------------------------------------------------------
# 3. Generate the theoretical power curve
# ------------------------------------------------------------

mu_grid <- seq(0, 4, length.out = 500)

power_grid <- pnorm(
  -z_critical + mu_grid / sigma
)

jpeg(
  filename = "power-curve-taustar-mult.jpg",
  width = 1800,
  height = 1200,
  res = 200,
  quality = 100
)

plot(
  mu_grid,
  power_grid,
  type = "l",
  lwd = 2,
  xlim = c(0, 4),
  ylim = c(0, 1),
  xlab = expression(mu),
  ylab = "Power",
  main = expression(
    "Theoretical power curve for " ~ hat(tau)^"*"
  )
)

abline(h = alpha, lty = 2)
points(
  mu_values,
  theoretical_power,
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

# ------------------------------------------------------------
# 4. Compare theoretical and simulated power
# ------------------------------------------------------------

jpeg(
  filename = "power-curve-taustar-simulation.jpg",
  width = 1800,
  height = 1200,
  res = 200,
  quality = 100
)

plot(
  mu_values,
  theoretical_power,
  type = "b",
  pch = 16,
  lty = 1,
  lwd = 2,
  xlim = c(0, 4),
  ylim = c(0, 1),
  xlab = expression(mu),
  ylab = "Power",
  main = "Theoretical and Monte Carlo Power"
)

lines(
  mu_values,
  simulated_power,
  type = "b",
  pch = 17,
  lty = 2,
  lwd = 2
)

abline(h = alpha, lty = 3)

legend(
  "bottomright",
  legend = c(
    "Theoretical power",
    "Monte Carlo power",
    "Significance level"
  ),
  pch = c(16, 17, NA),
  lty = c(1, 2, 3),
  lwd = c(2, 2, 1),
  bty = "n"
)

dev.off()