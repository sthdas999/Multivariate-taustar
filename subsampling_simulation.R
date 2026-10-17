###############################################################################
# SUBSAMPLING-BASED SIMULATION FOR A MULTIVARIATE RANK/SIGN DEPENDENCE STATISTIC
#
# Purpose:
#   1. Generate independent, linear, quadratic, and sinusoidal dependence data.
#   2. Compute the manuscript's custom multivariate fourth-order sign statistic.
#   3. Approximate the statistic by Monte Carlo sampling of quadruples, avoiding
#      the O(n^4) cost of enumerating every quadruple.
#   4. Construct subsampling distributions and root-n subsampling confidence
#      intervals.
#   5. Summarize empirical size/power, coverage, interval width, and quantiles.
#   6. Save CSV results and PNG plots.
#
# IMPORTANT METHODOLOGICAL NOTES:
#   * The kernel below reproduces the custom sign-product construction in the
#     manuscript's Parkinson/R-code section. It should NOT be called the
#     canonical Bergsma-Dassios tau* statistic until its kernel has been
#     mathematically verified against the definition/theorems in the paper.
#   * Random quadruple sampling approximates the all-quadruple average. It is
#     much faster, but adds Monte Carlo approximation error. Increase N_QUADS
#     and repeat with different seeds for a stability check.
#   * The confidence interval below uses sqrt(n) scaling, as in the manuscript's
#     subsampling section. Its validity depends on a suitable root-n limit.
#     Degenerate fourth-order U-statistics under independence may not have this
#     limit; check the theorem and kernel degeneracy before interpreting the
#     intervals as formally valid.
#   * Coverage for alternatives is computed against a large-sample Monte Carlo
#     reference for the same custom kernel. This is a numerical reference, not
#     an analytic population value. Under the null, the reference is set to 0.
###############################################################################

## ----------------------------- USER SETTINGS -------------------------------

set.seed(20261017)

OUTPUT_DIR <- "subsampling_simulation_results"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Quick test settings. Increase for a final study.
N_VALUES       <- c(60, 100, 150)
P               <- 3
R_MC            <- 100       # Monte Carlo data-generation replications
M_SUBSAMPLES    <- 100       # subsamples per generated dataset
N_QUADS         <- 10000     # sampled quadruples per statistic evaluation
N_QUADS_REF     <- 200000    # sampled quadruples for reference parameter
N_REF           <- 4000      # independent reference sample size
ALPHA           <- 0.05
SUBSAMPLE_POWER <- 0.70      # b_n = floor(n^gamma)
NOISE_SD        <- 0.70
LINEAR_RHO      <- 0.60
SINE_FREQUENCY  <- 1.50

# Use "quick" first. For final results, consider R_MC >= 500, M_SUBSAMPLES
# >= 500, and a sensitivity analysis over N_QUADS and SUBSAMPLE_POWER.
SCENARIOS <- c("independence", "linear", "quadratic", "sine")

## --------------------------- DATA GENERATION -------------------------------

generate_data <- function(n, p = 3, scenario,
                          noise_sd = 0.70, rho = 0.60,
                          sine_frequency = 1.50) {
  stopifnot(n >= 4, p >= 1)

  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", seq_len(p))
  eps <- rnorm(n)

  if (scenario == "independence") {
    # X and Y are generated independently.
    Y <- rnorm(n)
  } else if (scenario == "linear") {
    # Unit-variance linear signal when X1 and eps are standard normal.
    Y <- rho * X[, 1] + sqrt(1 - rho^2) * eps
  } else if (scenario == "quadratic") {
    # Nonlinear dependence with zero population Pearson correlation
    # between X1 and X1^2 + independent noise under symmetric X1.
    Y <- X[, 1]^2 + noise_sd * eps
  } else if (scenario == "sine") {
    # Oscillatory, non-monotone dependence.
    Y <- sin(sine_frequency * X[, 1]) + noise_sd * eps
  } else {
    stop("Unknown scenario: ", scenario)
  }

  list(X = X, Y = as.numeric(Y))
}

## ----------------------- CUSTOM SIGN-BASED STATISTIC -----------------------

# For a vector-valued row pair, this is the product of componentwise signs.
# This function is used in the custom kernel reproduced from the manuscript.
sign_product_rows <- function(X, a, b) {
  if (ncol(X) == 1L) return(sign(X[a, 1] - X[b, 1]))
  row_prod <- rep(1, length(a))
  for (j in seq_len(ncol(X))) {
    row_prod <- row_prod * sign(X[cbind(a, j)] - X[cbind(b, j)])
  }
  row_prod
}

# Monte Carlo approximation to the manuscript's custom fourth-order statistic.
# Each sampled quadruple consists of four distinct row indices. Their order
# follows increasing row index, matching the manuscript's combn(n, 4) pattern.
#
# The exact version would average this kernel over choose(n, 4) quadruples.
# Sampling N_QUADS quadruples is an approximation to that average.
tau_star_mv_mc <- function(X, Y, n_quads = 10000) {
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  n <- nrow(X)

  if (length(Y) != n) stop("nrow(X) must equal length(Y).")
  if (n < 4L) stop("At least four observations are required.")

  q <- min(as.integer(n_quads), choose(n, 4))
  if (q <= 0) stop("n_quads must be positive.")

  # If all combinations are small enough, enumerate them exactly.
  # Otherwise sample quadruples without replacement where feasible.
  total_quad <- choose(n, 4)
  if (total_quad <= q && total_quad < 2e6) {
    cmb <- utils::combn(n, 4)
    i <- cmb[1, ]; j <- cmb[2, ]; k <- cmb[3, ]; ell <- cmb[4, ]
  } else {
    # Generate distinct indices within each row, then sort each quadruple.
    # Duplicate quadruples across rows are allowed; this is Monte Carlo
    # sampling from the set of all quadruples.
    ind <- matrix(NA_integer_, nrow = q, ncol = 4)
    for (r in seq_len(q)) ind[r, ] <- sort(sample.int(n, 4L, replace = FALSE))
    i <- ind[, 1]; j <- ind[, 2]; k <- ind[, 3]; ell <- ind[, 4]
  }

  # Kernel as coded in the manuscript's custom Parkinson example:
  # sign_product(X_i, X_j) * sign_product(X_k, X_ell)
  # multiplied by sign(Y_i-Y_j) * sign(Y_k-Y_ell).
  sx1 <- sign_product_rows(X, i, j)
  sx2 <- sign_product_rows(X, k, ell)
  sy1 <- sign(Y[i] - Y[j])
  sy2 <- sign(Y[k] - Y[ell])

  mean(sx1 * sx2 * sy1 * sy2)
}

## ------------------------- SUBSAMPLING HELPERS ------------------------------

get_subsample_size <- function(n, gamma = 0.70) {
  b <- floor(n^gamma)
  b <- max(4L, min(as.integer(b), n - 1L))
  b
}

subsampling_distribution <- function(X, Y, tau_full, b, M = 100,
                                     n_quads = 10000) {
  n <- nrow(X)
  if (b >= n) stop("Subsample size b must be smaller than n.")

  vals <- numeric(M)
  for (m in seq_len(M)) {
    idx <- sample.int(n, size = b, replace = FALSE)
    tau_b <- tau_star_mv_mc(X[idx, , drop = FALSE], Y[idx],
                            n_quads = n_quads)
    vals[m] <- sqrt(b) * (tau_b - tau_full)
  }
  vals
}

# Basic two-sided subsampling CI obtained by quantile inversion:
# [tau_hat - q_(1-alpha/2)/sqrt(n), tau_hat - q_(alpha/2)/sqrt(n)].
# Use only if the relevant root-n limit is justified for the statistic.
subsampling_ci <- function(tau_hat, sub_dist, n, alpha = 0.05) {
  qq <- quantile(sub_dist,
                 probs = c(alpha / 2, 1 - alpha / 2),
                 names = FALSE, type = 8)
  lower <- tau_hat - qq[2] / sqrt(n)
  upper <- tau_hat - qq[1] / sqrt(n)
  c(lower = lower, upper = upper)
}

## --------------------- POPULATION-REFERENCE APPROXIMATION ------------------

# For alternatives, estimate the population target of the custom kernel using
# an independent, large sample and many sampled quadruples. Under independence,
# the target is set to zero for the stated construction.
estimate_reference_tau <- function(scenario, p, n_ref, n_quads_ref,
                                   noise_sd, rho, sine_frequency,
                                   seed = 1234) {
  if (scenario == "independence") return(0)

  set.seed(seed)
  ref <- generate_data(n_ref, p, scenario, noise_sd, rho, sine_frequency)
  tau_star_mv_mc(ref$X, ref$Y, n_quads = n_quads_ref)
}

## ------------------------- ONE DATASET REPLICATION --------------------------

run_one_replication <- function(n, scenario, p, M, n_quads, gamma, alpha,
                                tau_reference, noise_sd, rho,
                                sine_frequency) {
  dat <- generate_data(n, p, scenario, noise_sd, rho, sine_frequency)
  tau_hat <- tau_star_mv_mc(dat$X, dat$Y, n_quads = n_quads)
  b <- get_subsample_size(n, gamma)

  sub_dist <- subsampling_distribution(
    X = dat$X, Y = dat$Y, tau_full = tau_hat, b = b,
    M = M, n_quads = n_quads
  )

  ci <- subsampling_ci(tau_hat, sub_dist, n, alpha)
  reject_zero <- (ci["lower"] > 0 || ci["upper"] < 0)
  covered <- (ci["lower"] <= tau_reference &&
              tau_reference <= ci["upper"])

  data.frame(
    scenario = scenario,
    n = n,
    p = p,
    b = b,
    tau_hat = tau_hat,
    tau_reference = tau_reference,
    ci_lower = unname(ci["lower"]),
    ci_upper = unname(ci["upper"]),
    ci_width = unname(ci["upper"] - ci["lower"]),
    covers_reference = as.integer(covered),
    reject_zero = as.integer(reject_zero),
    sub_q025 = unname(quantile(sub_dist, alpha / 2, names = FALSE)),
    sub_median = unname(quantile(sub_dist, 0.50, names = FALSE)),
    sub_q975 = unname(quantile(sub_dist, 1 - alpha / 2, names = FALSE))
  )
}

## ------------------------------ RUN STUDY -----------------------------------

all_results <- list()
result_counter <- 1L

for (scenario in SCENARIOS) {
  message("\nScenario: ", scenario)

  tau_ref <- estimate_reference_tau(
    scenario = scenario,
    p = P,
    n_ref = N_REF,
    n_quads_ref = N_QUADS_REF,
    noise_sd = NOISE_SD,
    rho = LINEAR_RHO,
    sine_frequency = SINE_FREQUENCY,
    seed = 7000 + match(scenario, SCENARIOS)
  )

  message("Reference value for custom kernel: ",
          format(tau_ref, digits = 5))

  for (n in N_VALUES) {
    message("  n = ", n, "; b = ", get_subsample_size(n, SUBSAMPLE_POWER))

    rep_results <- vector("list", R_MC)

    for (r in seq_len(R_MC)) {
      # Reproducible but different random stream for each replication.
      set.seed(100000 + 1000 * match(scenario, SCENARIOS) + 10 * n + r)

      rep_results[[r]] <- run_one_replication(
        n = n,
        scenario = scenario,
        p = P,
        M = M_SUBSAMPLES,
        n_quads = N_QUADS,
        gamma = SUBSAMPLE_POWER,
        alpha = ALPHA,
        tau_reference = tau_ref,
        noise_sd = NOISE_SD,
        rho = LINEAR_RHO,
        sine_frequency = SINE_FREQUENCY
      )

      if (r %% 10 == 0) {
        message("    completed replication ", r, " / ", R_MC)
      }
    }

    all_results[[result_counter]] <- do.call(rbind, rep_results)
    result_counter <- result_counter + 1L
  }
}

replication_results <- do.call(rbind, all_results)
row.names(replication_results) <- NULL

write.csv(
  replication_results,
  file = file.path(OUTPUT_DIR, "subsampling_replication_results.csv"),
  row.names = FALSE
)

## ---------------------------- SUMMARY TABLES -------------------------------

# The reported "rejection probability" is the fraction of intervals excluding
# zero. Under the independence null, this estimates empirical size. Under
# alternatives it is an interval-based rejection rate/power proxy.
summarise_group <- function(d) {
  data.frame(
    scenario = d$scenario[1],
    n = d$n[1],
    p = d$p[1],
    b = d$b[1],
    mean_tau_hat = mean(d$tau_hat, na.rm = TRUE),
    sd_tau_hat = sd(d$tau_hat, na.rm = TRUE),
    empirical_rejection = mean(d$reject_zero, na.rm = TRUE),
    rejection_mcse = sqrt(
      mean(d$reject_zero, na.rm = TRUE) *
        (1 - mean(d$reject_zero, na.rm = TRUE)) / nrow(d)
    ),
    empirical_coverage = mean(d$covers_reference, na.rm = TRUE),
    mean_ci_width = mean(d$ci_width, na.rm = TRUE),
    sd_ci_width = sd(d$ci_width, na.rm = TRUE),
    reference_tau = d$tau_reference[1]
  )
}

summary_list <- split(
  replication_results,
  interaction(replication_results$scenario,
              replication_results$n, drop = TRUE)
)
summary_results <- do.call(rbind, lapply(summary_list, summarise_group))
row.names(summary_results) <- NULL

write.csv(
  summary_results,
  file = file.path(OUTPUT_DIR, "subsampling_summary_results.csv"),
  row.names = FALSE
)

cat("\n================ SUBSAMPLING SUMMARY ================\n")
print(summary_results, digits = 4, row.names = FALSE)

## ------------------------------ DIAGNOSTIC PLOTS ---------------------------

# 1. Empirical rejection rate by sample size and scenario.
png(file.path(OUTPUT_DIR, "empirical_rejection.png"),
    width = 1200, height = 800, res = 140)
scens <- unique(summary_results$scenario)
plot(NA, xlim = range(N_VALUES), ylim = c(0, 1),
     xlab = "Sample size n", ylab = "Empirical rejection probability",
     main = "Subsampling interval test: empirical rejection")
for (s in scens) {
  ds <- summary_results[summary_results$scenario == s, ]
  lines(ds$n, ds$empirical_rejection, type = "b", pch = 19)
}
abline(h = ALPHA, lty = 2)
legend("topleft", legend = c(scens, paste0("Nominal alpha = ", ALPHA)),
       lty = c(rep(1, length(scens)), 2),
       pch = c(rep(19, length(scens)), NA), bty = "n")
dev.off()

# 2. Empirical coverage by sample size and scenario.
png(file.path(OUTPUT_DIR, "empirical_coverage.png"),
    width = 1200, height = 800, res = 140)
plot(NA, xlim = range(N_VALUES), ylim = c(0, 1),
     xlab = "Sample size n", ylab = "Empirical coverage probability",
     main = "Subsampling confidence interval coverage")
for (s in scens) {
  ds <- summary_results[summary_results$scenario == s, ]
  lines(ds$n, ds$empirical_coverage, type = "b", pch = 19)
}
abline(h = 1 - ALPHA, lty = 2)
legend("bottomright", legend = c(scens, paste0("Nominal coverage = ", 1 - ALPHA)),
       lty = c(rep(1, length(scens)), 2),
       pch = c(rep(19, length(scens)), NA), bty = "n")
dev.off()

# 3. Average interval width by sample size and scenario.
png(file.path(OUTPUT_DIR, "mean_interval_width.png"),
    width = 1200, height = 800, res = 140)
plot(NA, xlim = range(N_VALUES),
     ylim = range(summary_results$mean_ci_width, finite = TRUE),
     xlab = "Sample size n", ylab = "Mean interval width",
     main = "Average width of subsampling confidence intervals")
for (s in scens) {
  ds <- summary_results[summary_results$scenario == s, ]
  lines(ds$n, ds$mean_ci_width, type = "b", pch = 19)
}
legend("topright", legend = scens, lty = 1, pch = 19, bty = "n")
dev.off()

# 4. Distribution of full-sample estimates by scenario and n.
png(file.path(OUTPUT_DIR, "tau_estimate_boxplots.png"),
    width = 1400, height = 900, res = 140)
boxplot(tau_hat ~ interaction(scenario, n),
        data = replication_results, las = 2,
        xlab = "Scenario and sample size",
        ylab = expression(hat(tau)^"*"),
        main = "Distribution of full-sample dependence estimates")
dev.off()

cat("\nFiles written to: ", normalizePath(OUTPUT_DIR), "\n", sep = "")
cat(" - subsampling_replication_results.csv\n")
cat(" - subsampling_summary_results.csv\n")
cat(" - empirical_rejection.png\n")
cat(" - empirical_coverage.png\n")
cat(" - mean_interval_width.png\n")
cat(" - tau_estimate_boxplots.png\n")

###############################################################################
# OPTIONAL SENSITIVITY ANALYSIS
#
# Repeat the study with:
#   SUBSAMPLE_POWER <- 0.60 and 0.80
#   N_QUADS <- 5000, 10000, 25000
#   M_SUBSAMPLES <- 200 or 500
#   R_MC <- 500 or more
#
# Do not compare results across settings unless the random-quadruple
# approximation error has been checked. For final publication, validate the
# kernel mathematically and compare the Monte Carlo approximation with exact
# enumeration for small n (e.g., n <= 30).
###############################################################################
