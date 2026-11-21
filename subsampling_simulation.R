###############################################################################
# MULTIVARIATE FOURTH-ORDER SIGN-CONCORDANCE STATISTIC
# Revised implementation based on Sections 3, 4 and 5 of the manuscript
#
# Main features:
#  1. Uses the concordance comparison sign{(x1-x3)(x2-x4)} in Section 3.1.
#  2. Symmetrizes the fourth-order kernel over its three distinct pairings.
#  3. Computes the fourth-order U-statistic by exact enumeration when feasible,
#     and Monte Carlo quadruple sampling otherwise.
#  4. Implements the centred, sqrt(b_n)-scaled subsampling distribution
#     described in Section 5.1.
#  5. Constructs the manuscript's two-sided confidence interval.
#  6. Implements the upper-tail subsampling test described in Section 5.1.
#  7. Reports interval-based rejection separately from the upper-tail test.
#  8. Estimates scenario-specific reference values by independent simulation.
#
# THEORETICAL CAUTION:
# The manuscript's root-n Gaussian limit and second-order degeneracy assertion
# require a separate mathematical verification. The sqrt(n) normalization
# below reproduces the proposed procedure; it does not prove its validity.
###############################################################################

set.seed(20261121)

## 1. USER SETTINGS -----------------------------------------------------------

OUTPUT_DIR <- "revised_multivariate_taustar_results"
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

N_VALUES        <- c(60, 100, 150)
P               <- 3
Q               <- 1

R_MC            <- 100
M_SUBSAMPLES    <- 300

N_QUADS_FULL    <- 100000
N_QUADS_REF     <- 200000
N_REF           <- 4000

ALPHA           <- 0.05
SUBSAMPLE_POWER <- 0.70

NOISE_SD        <- 0.70
LINEAR_RHO      <- 0.60
SINE_FREQUENCY  <- 1.50

# Exact enumeration is used whenever choose(n,4) does not exceed this limit.
EXACT_LIMIT     <- 2000000

SCENARIOS <- c("independence", "linear", "quadratic", "sine")


## 2. DATA GENERATION ---------------------------------------------------------

generate_data <- function(n, p = 3, q = 1,
                          scenario = c("independence", "linear",
                                       "quadratic", "sine"),
                          noise_sd = 0.70,
                          rho = 0.60,
                          sine_frequency = 1.50) {

  scenario <- match.arg(scenario)

  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", seq_len(p))

  Y <- matrix(rnorm(n * q), nrow = n, ncol = q)
  colnames(Y) <- paste0("Y", seq_len(q))

  eps <- rnorm(n)

  if (scenario == "independence") {

    # X and Y are generated independently.

  } else if (scenario == "linear") {

    Y[, 1] <- rho * X[, 1] +
      sqrt(1 - rho^2) * eps

  } else if (scenario == "quadratic") {

    Y[, 1] <- X[, 1]^2 + noise_sd * eps

  } else if (scenario == "sine") {

    Y[, 1] <- sin(sine_frequency * X[, 1]) +
      noise_sd * eps
  }

  list(X = X, Y = Y)
}


## 3. MULTIVARIATE CONCORDANCE COMPONENT -------------------------------------

# For a particular pairing (r,s) and (t,u), calculate
#
#   product_j sign{(X_rj - X_sj)(X_tj - X_uj)}.
#
# This is the coordinatewise extension of the comparison in Section 3.1.

concordance_component <- function(A, i1, i2, i3, i4) {

  A <- as.matrix(A)

  ans <- rep(1, length(i1))

  for (j in seq_len(ncol(A))) {

    ans <- ans *
      sign((A[i1, j] - A[i2, j]) *
           (A[i3, j] - A[i4, j]))
  }

  ans
}


## 4. SYMMETRIC FOURTH-ORDER KERNEL ------------------------------------------

# The displayed comparison in Section 3.1 uses one pairing:
# (1,3) with (2,4).
#
# A fourth-order U-statistic requires a symmetric kernel. The code therefore
# uses the average over the three distinct pairings of four observations.
#
# This is the symmetric version of the product construction. It should be
# checked against the precise intended population definition before being
# identified with the canonical Bergsma-Dassios tau-star kernel.

kernel_tau4 <- function(X, Y, ind) {

  X <- as.matrix(X)
  Y <- as.matrix(Y)
  ind <- as.matrix(ind)

  i <- ind[, 1]
  j <- ind[, 2]
  k <- ind[, 3]
  l <- ind[, 4]

  # Pairing 1: (1,2) with (3,4)
  ax1 <- concordance_component(X, i, j, k, l)
  ay1 <- concordance_component(Y, i, j, k, l)

  # Pairing 2: (1,3) with (2,4)
  ax2 <- concordance_component(X, i, k, j, l)
  ay2 <- concordance_component(Y, i, k, j, l)

  # Pairing 3: (1,4) with (2,3)
  ax3 <- concordance_component(X, i, l, j, k)
  ay3 <- concordance_component(Y, i, l, j, k)

  h1 <- ax1 * ay1
  h2 <- ax2 * ay2
  h3 <- ax3 * ay3

  (h1 + h2 + h3) / 3
}


## 5. GENERATE RANDOM DISTINCT QUADRUPLES ------------------------------------

sample_quadruples <- function(n, B) {

  out <- matrix(NA_integer_, nrow = B, ncol = 4)

  for (r in seq_len(B)) {
    out[r, ] <- sort(sample.int(n, size = 4, replace = FALSE))
  }

  out
}


## 6. FOURTH-ORDER U-STATISTIC -----------------------------------------------

tau_star_mv <- function(X, Y,
                        n_quads = 100000,
                        exact_limit = 2000000) {

  X <- as.matrix(X)
  Y <- as.matrix(Y)

  n <- nrow(X)

  if (nrow(Y) != n) {
    stop("X and Y must have the same number of observations.")
  }

  if (n < 4) {
    stop("At least four observations are required.")
  }

  if (n_quads < 1) {
    stop("n_quads must be positive.")
  }

  total <- choose(n, 4)

  if (total <= exact_limit) {

    # Exact U-statistic.
    ind <- utils::combn(n, 4)
    ind <- t(ind)

    values <- kernel_tau4(X, Y, ind)

    return(mean(values))

  } else {

    # Monte Carlo approximation to the complete U-statistic.
    B <- min(as.integer(n_quads), as.integer(total))

    ind <- sample_quadruples(n, B)

    values <- kernel_tau4(X, Y, ind)

    return(mean(values))
  }
}


## 7. SUBSAMPLE SIZE ----------------------------------------------------------

get_subsample_size <- function(n, gamma = 0.70) {

  b <- floor(n^gamma)

  b <- max(4L, min(as.integer(b), n - 1L))

  b
}


## 8. SUBSAMPLING DISTRIBUTION ------------------------------------------------

# Section 5.1:
#
# T_m = sqrt(b_n) * (tau_hat_b(m) - tau_hat_n).
#
# The subsample estimates are calculated exactly whenever feasible.

subsampling_distribution <- function(X, Y, tau_full,
                                     b,
                                     M = 300,
                                     n_quads = 100000,
                                     exact_limit = 2000000) {

  n <- nrow(X)

  if (b >= n) {
    stop("Subsample size must be smaller than the original sample size.")
  }

  vals <- numeric(M)

  for (m in seq_len(M)) {

    idx <- sample.int(n, size = b, replace = FALSE)

    tau_b <- tau_star_mv(
      X[idx, , drop = FALSE],
      Y[idx, , drop = FALSE],
      n_quads = n_quads,
      exact_limit = exact_limit
    )

    vals[m] <- sqrt(b) * (tau_b - tau_full)
  }

  vals
}


## 9. SUBSAMPLING CONFIDENCE INTERVAL ----------------------------------------

# The confidence interval follows the quantile inversion in Section 5.1:
#
# [tau_hat - q_(1-alpha/2)/sqrt(n),
#  tau_hat - q_(alpha/2)/sqrt(n)].

subsampling_ci <- function(tau_hat, sub_dist,
                           n, alpha = 0.05) {

  qq <- quantile(
    sub_dist,
    probs = c(alpha / 2, 1 - alpha / 2),
    names = FALSE,
    type = 8
  )

  lower <- tau_hat - qq[2] / sqrt(n)
  upper <- tau_hat - qq[1] / sqrt(n)

  c(lower = lower, upper = upper)
}


## 10. UPPER-TAIL SUBSAMPLING TEST -------------------------------------------

# Section 5.1 proposes rejecting H0 when
#
# sqrt(n) * tau_hat > q_(1-alpha),
#
# where q_(1-alpha) is the upper quantile of the centred subsampling
# distribution.

subsampling_test <- function(tau_hat, sub_dist,
                             n, alpha = 0.05) {

  observed <- sqrt(n) * tau_hat

  critical_value <- as.numeric(
    quantile(sub_dist, 1 - alpha, names = FALSE, type = 8)
  )

  p_value <- (1 + sum(sub_dist >= observed)) /
    (length(sub_dist) + 1)

  list(
    statistic = observed,
    critical_value = critical_value,
    p_value = p_value,
    reject = as.integer(observed > critical_value)
  )
}


## 11. REFERENCE VALUE --------------------------------------------------------

# Under independence, the population target is zero according to the
# manuscript's population formulation.
#
# Under alternatives, use an independent large-sample Monte Carlo estimate.
# This is a numerical reference, not an analytic population value.

estimate_reference_tau <- function(
    scenario,
    p,
    q,
    n_ref,
    n_quads_ref,
    noise_sd,
    rho,
    sine_frequency,
    seed = 1234) {

  if (scenario == "independence") {
    return(0)
  }

  set.seed(seed)

  ref <- generate_data(
    n = n_ref,
    p = p,
    q = q,
    scenario = scenario,
    noise_sd = noise_sd,
    rho = rho,
    sine_frequency = sine_frequency
  )

  tau_star_mv(
    ref$X,
    ref$Y,
    n_quads = n_quads_ref,
    exact_limit = 0
  )
}


## 12. ONE MONTE CARLO REPLICATION -------------------------------------------

run_one_replication <- function(
    n,
    scenario,
    p,
    q,
    M,
    n_quads,
    gamma,
    alpha,
    tau_reference,
    noise_sd,
    rho,
    sine_frequency,
    exact_limit) {

  dat <- generate_data(
    n = n,
    p = p,
    q = q,
    scenario = scenario,
    noise_sd = noise_sd,
    rho = rho,
    sine_frequency = sine_frequency
  )

  tau_hat <- tau_star_mv(
    dat$X,
    dat$Y,
    n_quads = n_quads,
    exact_limit = exact_limit
  )

  b <- get_subsample_size(n, gamma)

  sub_dist <- subsampling_distribution(
    X = dat$X,
    Y = dat$Y,
    tau_full = tau_hat,
    b = b,
    M = M,
    n_quads = n_quads,
    exact_limit = exact_limit
  )

  ci <- subsampling_ci(
    tau_hat = tau_hat,
    sub_dist = sub_dist,
    n = n,
    alpha = alpha
  )

  test <- subsampling_test(
    tau_hat = tau_hat,
    sub_dist = sub_dist,
    n = n,
    alpha = alpha
  )

  # Interval-based rejection, as used in the original simulation section.
  reject_interval <- as.integer(
    ci["lower"] > 0 || ci["upper"] < 0
  )

  covered <- as.integer(
    ci["lower"] <= tau_reference &&
      tau_reference <= ci["upper"]
  )

  data.frame(
    scenario = scenario,
    n = n,
    p = p,
    q = q,
    b = b,

    tau_hat = tau_hat,
    tau_reference = tau_reference,

    ci_lower = unname(ci["lower"]),
    ci_upper = unname(ci["upper"]),
    ci_width = unname(ci["upper"] - ci["lower"]),

    covers_reference = covered,

    # Keep the two decision rules separate.
    reject_interval = reject_interval,
    reject_upper_tail = test$reject,

    upper_tail_statistic = test$statistic,
    upper_tail_critical = test$critical_value,
    upper_tail_pvalue = test$p_value,

    sub_q025 = unname(
      quantile(sub_dist, alpha / 2, names = FALSE)
    ),

    sub_median = unname(
      quantile(sub_dist, 0.50, names = FALSE)
    ),

    sub_q975 = unname(
      quantile(sub_dist, 1 - alpha / 2, names = FALSE)
    )
  )
}


## 13. RUN THE SIMULATION -----------------------------------------------------

all_results <- list()
counter <- 1L

for (scenario in SCENARIOS) {

  message("\nScenario: ", scenario)

  tau_ref <- estimate_reference_tau(
    scenario = scenario,
    p = P,
    q = Q,
    n_ref = N_REF,
    n_quads_ref = N_QUADS_REF,
    noise_sd = NOISE_SD,
    rho = LINEAR_RHO,
    sine_frequency = SINE_FREQUENCY,
    seed = 7000 + match(scenario, SCENARIOS)
  )

  message(
    "Reference value = ",
    format(tau_ref, digits = 6)
  )

  for (n in N_VALUES) {

    b <- get_subsample_size(n, SUBSAMPLE_POWER)

    message("  n = ", n, "; b = ", b)

    rep_results <- vector("list", R_MC)

    for (r in seq_len(R_MC)) {

      set.seed(
        100000 +
          10000 * match(scenario, SCENARIOS) +
          100 * n + r
      )

      rep_results[[r]] <- run_one_replication(
        n = n,
        scenario = scenario,
        p = P,
        q = Q,
        M = M_SUBSAMPLES,
        n_quads = N_QUADS_FULL,
        gamma = SUBSAMPLE_POWER,
        alpha = ALPHA,
        tau_reference = tau_ref,
        noise_sd = NOISE_SD,
        rho = LINEAR_RHO,
        sine_frequency = SINE_FREQUENCY,
        exact_limit = EXACT_LIMIT
      )

      if (r %% 10 == 0) {
        message("    Replication ", r, " / ", R_MC)
      }
    }

    all_results[[counter]] <- do.call(rbind, rep_results)

    counter <- counter + 1L
  }
}

replication_results <- do.call(rbind, all_results)
row.names(replication_results) <- NULL

write.csv(
  replication_results,
  file = file.path(
    OUTPUT_DIR,
    "subsampling_replication_results.csv"
  ),
  row.names = FALSE
)


## 14. SUMMARY TABLE ----------------------------------------------------------

summarise_group <- function(d) {

  rejection_interval <- mean(d$reject_interval)
  rejection_upper <- mean(d$reject_upper_tail)

  data.frame(
    scenario = d$scenario[1],
    n = d$n[1],
    p = d$p[1],
    q = d$q[1],
    b = d$b[1],

    mean_tau_hat = mean(d$tau_hat),
    sd_tau_hat = sd(d$tau_hat),

    interval_rejection = rejection_interval,

    interval_rejection_mcse = sqrt(
      rejection_interval * (1 - rejection_interval) / nrow(d)
    ),

    upper_tail_rejection = rejection_upper,

    upper_tail_rejection_mcse = sqrt(
      rejection_upper * (1 - rejection_upper) / nrow(d)
    ),

    empirical_coverage = mean(d$covers_reference),

    coverage_mcse = sqrt(
      mean(d$covers_reference) *
        (1 - mean(d$covers_reference)) / nrow(d)
    ),

    mean_ci_width = mean(d$ci_width),
    sd_ci_width = sd(d$ci_width),

    mean_upper_tail_pvalue = mean(d$upper_tail_pvalue),

    reference_tau = d$tau_reference[1]
  )
}

summary_list <- split(
  replication_results,
  interaction(
    replication_results$scenario,
    replication_results$n,
    drop = TRUE
  )
)

summary_results <- do.call(
  rbind,
  lapply(summary_list, summarise_group)
)

row.names(summary_results) <- NULL

write.csv(
  summary_results,
  file = file.path(
    OUTPUT_DIR,
    "subsampling_summary_results.csv"
  ),
  row.names = FALSE
)

cat("\n================ SUMMARY RESULTS ================\n")
print(summary_results, digits = 4, row.names = FALSE)


## 15. DIAGNOSTIC PLOTS -------------------------------------------------------

scens <- unique(summary_results$scenario)

# Plot 1: Interval-based rejection rates.

png(
  file.path(OUTPUT_DIR, "interval_rejection_rates.png"),
  width = 1200, height = 800, res = 140
)

plot(
  NA,
  xlim = range(N_VALUES),
  ylim = c(0, 1),
  xlab = "Sample size n",
  ylab = "Interval-based rejection probability",
  main = "Subsampling confidence-interval rejection rates"
)

for (s in scens) {
  ds <- summary_results[summary_results$scenario == s, ]
  lines(ds$n, ds$interval_rejection, type = "b", pch = 19)
}

abline(h = ALPHA, lty = 2)

legend(
  "topleft",
  legend = c(scens, paste0("Nominal level = ", ALPHA)),
  lty = c(rep(1, length(scens)), 2),
  pch = c(rep(19, length(scens)), NA),
  bty = "n"
)

dev.off()


# Plot 2: Upper-tail test rejection rates.

png(
  file.path(OUTPUT_DIR, "upper_tail_rejection_rates.png"),
  width = 1200, height = 800, res = 140
)

plot(
  NA,
  xlim = range(N_VALUES),
  ylim = c(0, 1),
  xlab = "Sample size n",
  ylab = "Upper-tail rejection probability",
  main = "Upper-tail subsampling test"
)

for (s in scens) {
  ds <- summary_results[summary_results$scenario == s, ]
  lines(ds$n, ds$upper_tail_rejection, type = "b", pch = 19)
}

abline(h = ALPHA, lty = 2)

legend(
  "topleft",
  legend = c(scens, paste0("Nominal level = ", ALPHA)),
  lty = c(rep(1, length(scens)), 2),
  pch = c(rep(19, length(scens)), NA),
  bty = "n"
)

dev.off()


# Plot 3: Coverage probabilities.

png(
  file.path(OUTPUT_DIR, "empirical_coverage.png"),
  width = 1200, height = 800, res = 140
)

plot(
  NA,
  xlim = range(N_VALUES),
  ylim = c(0, 1),
  xlab = "Sample size n",
  ylab = "Empirical coverage probability",
  main = "Subsampling confidence-interval coverage"
)

for (s in scens) {
  ds <- summary_results[summary_results$scenario == s, ]
  lines(ds$n, ds$empirical_coverage, type = "b", pch = 19)
}

abline(h = 1 - ALPHA, lty = 2)

legend(
  "bottomright",
  legend = c(scens, paste0("Nominal coverage = ", 1 - ALPHA)),
  lty = c(rep(1, length(scens)), 2),
  pch = c(rep(19, length(scens)), NA),
  bty = "n"
)

dev.off()


# Plot 4: Mean interval width.

png(
  file.path(OUTPUT_DIR, "mean_interval_width.png"),
  width = 1200, height = 800, res = 140
)

plot(
  NA,
  xlim = range(N_VALUES),
  ylim = range(summary_results$mean_ci_width, finite = TRUE),
  xlab = "Sample size n",
  ylab = "Mean interval width",
  main = "Average subsampling confidence-interval width"
)

for (s in scens) {
  ds <- summary_results[summary_results$scenario == s, ]
  lines(ds$n, ds$mean_ci_width, type = "b", pch = 19)
}

legend(
  "topright",
  legend = scens,
  lty = 1,
  pch = 19,
  bty = "n"
)

dev.off()


# Plot 5: Distribution of estimates.

png(
  file.path(OUTPUT_DIR, "tau_estimate_boxplots.png"),
  width = 1400, height = 900, res = 140
)

boxplot(
  tau_hat ~ interaction(scenario, n),
  data = replication_results,
  las = 2,
  xlab = "Scenario and sample size",
  ylab = expression(hat(tau)^"*"),
  main = "Distribution of full-sample estimates"
)

dev.off()


## 16. SAVE SETTINGS AND COMPLETION MESSAGE ----------------------------------

settings <- data.frame(
  parameter = c(
    "N_VALUES", "P", "Q", "R_MC", "M_SUBSAMPLES",
    "N_QUADS_FULL", "N_QUADS_REF", "N_REF",
    "ALPHA", "SUBSAMPLE_POWER", "EXACT_LIMIT"
  ),
  value = c(
    paste(N_VALUES, collapse = ","),
    P, Q, R_MC, M_SUBSAMPLES,
    N_QUADS_FULL, N_QUADS_REF, N_REF,
    ALPHA, SUBSAMPLE_POWER, EXACT_LIMIT
  )
)

write.csv(
  settings,
  file = file.path(OUTPUT_DIR, "simulation_settings.csv"),
  row.names = FALSE
)

cat("\nSimulation completed.\n")
cat("Results directory: ", normalizePath(OUTPUT_DIR), "\n", sep = "")

###############################################################################
# END OF PROGRAM
###############################################################################
