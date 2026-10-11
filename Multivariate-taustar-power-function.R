# ============================================================
# Numerical results for the independence test
# Section 3: Proposed Test Statistic
# ============================================================

set.seed(123)

sample_sizes <- c(50, 100, 200, 300)
B <- 1000
M <- 50000
alpha <- 0.05
rho <- 0.5

# ------------------------------------------------------------
# Multivariate concordance function defined in the manuscript
# ------------------------------------------------------------

a_quad <- function(A) {
  # A is an array with dimensions 4 x p x M
  p <- dim(A)[2]
  m <- dim(A)[3]
  ans <- rep(1, m)
  
  for (j in seq_len(p)) {
    z <- (A[1, j, ] - A[3, j, ]) *
      (A[2, j, ] - A[4, j, ])
    
    ans <- ans * sign(z)
  }
  
  ans
}

# ------------------------------------------------------------
# Monte Carlo approximation to the order-four U-statistic
# ------------------------------------------------------------

tau_star <- function(X, Y, M = 50000) {
  n <- nrow(X)
  
  ids <- matrix(
    sample.int(n, 4 * M, replace = TRUE),
    nrow = 4
  )
  
  # Require four distinct observations per quadruple
  bad <- apply(ids, 2, function(z) length(unique(z)) < 4)
  
  while (any(bad)) {
    ids[, bad] <- matrix(
      sample.int(n, 4 * sum(bad), replace = TRUE),
      nrow = 4
    )
    bad <- apply(ids, 2, function(z) length(unique(z)) < 4)
  }
  
  Xarr <- array(NA_real_, dim = c(4, ncol(X), M))
  Yarr <- array(NA_real_, dim = c(4, ncol(Y), M))
  
  for (k in 1:4) {
    Xarr[k, , ] <- t(X[ids[k, ], , drop = FALSE])
    Yarr[k, , ] <- t(Y[ids[k, ], , drop = FALSE])
  }
  
  mean(a_quad(Xarr) * a_quad(Yarr))
}

# ------------------------------------------------------------
# Monte Carlo estimate of the second-order projection variance
# sigma^2 = 4 Var{E[h(Z1,Z2,Z3,Z4) | Z1,Z2]}
# ------------------------------------------------------------

estimate_variance <- function(X, Y, Mvar = 20000) {
  n <- nrow(X)
  
  pair_ids <- matrix(
    sample.int(n, 2 * Mvar, replace = TRUE),
    nrow = 2
  )
  
  other_ids <- matrix(
    sample.int(n, 2 * Mvar, replace = TRUE),
    nrow = 2
  )
  
  h2 <- numeric(Mvar)
  
  for (k in seq_len(Mvar)) {
    ids <- c(pair_ids[, k], other_ids[, k])
    
    if (length(unique(ids)) < 4) {
      h2[k] <- NA_real_
      next
    }
    
    xarr <- array(NA_real_, dim = c(4, ncol(X), 1))
    yarr <- array(NA_real_, dim = c(4, ncol(Y), 1))
    
    for (j in 1:4) {
      xarr[j, , 1] <- X[ids[j], ]
      yarr[j, , 1] <- Y[ids[j], ]
    }
    
    h2[k] <- a_quad(xarr) * a_quad(yarr)
  }
  
  h2 <- h2[is.finite(h2)]
  
  if (length(h2) < 2) return(NA_real_)
  
  4 * var(h2)
}

# ------------------------------------------------------------
# Independence test: permutation and asymptotic procedures
# ------------------------------------------------------------

independence_test <- function(n, B = 1000, M = 50000) {
  
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  
  dat <- MASS::mvrnorm(
    n = n,
    mu = c(0, 0),
    Sigma = Sigma
  )
  
  X <- matrix(dat[, 1], ncol = 1)
  Y <- matrix(dat[, 2], ncol = 1)
  
  tau_obs <- tau_star(X, Y, M = M)
  
  sigma2_hat <- estimate_variance(X, Y)
  
  if (is.finite(sigma2_hat) && sigma2_hat > 0) {
    z_stat <- sqrt(n) * tau_obs / sqrt(sigma2_hat)
    p_asym <- 2 * pnorm(-abs(z_stat))
  } else {
    p_asym <- NA_real_
  }
  
  tau_perm <- numeric(B)
  
  for (b in seq_len(B)) {
    Y_perm <- Y[sample.int(n), , drop = FALSE]
    tau_perm[b] <- tau_star(X, Y_perm, M = M)
  }
  
  p_perm <- (1 + sum(abs(tau_perm) >= abs(tau_obs))) / (B + 1)
  
  data.frame(
    Sample_Size = n,
    Tau_Star = tau_obs,
    Sigma_Squared = sigma2_hat,
    Asymptotic_p = p_asym,
    Permutation_p = p_perm,
    Asymptotic_Decision = ifelse(
      is.na(p_asym), NA,
      ifelse(p_asym <= alpha, "Reject H0", "Do not reject H0")
    ),
    Permutation_Decision = ifelse(
      p_perm <= alpha, "Reject H0", "Do not reject H0"
    )
  )
}

# ------------------------------------------------------------
# Run the test for all sample sizes
# ------------------------------------------------------------

results <- do.call(
  rbind,
  lapply(
    sample_sizes,
    independence_test,
    B = B,
    M = M
  )
)

results[, 2:5] <- round(results[, 2:5], 4)

print(results, row.names = FALSE)