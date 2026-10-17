# ============================================================
# Numerical results for the independence test
# Section 3: Proposed Test Statistic
# Display generated results at every step
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, digits = 6)

set.seed(123)

sample_sizes <- c(50, 100, 200, 300)
B <- 1000
M <- 50000
Mvar <- 20000
alpha <- 0.05
rho <- 0.5

# Set TRUE to display detailed intermediate results
show_details <- TRUE

# ============================================================
# STEP 1: Multivariate concordance function
# ============================================================

a_quad <- function(A) {
  
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

cat("\n============================================\n")
cat("STEP 1: Concordance function defined\n")
cat("============================================\n")

# Example verification
A_example <- array(
  rnorm(4 * 1 * 5),
  dim = c(4, 1, 5)
)

cat("\nExample input array dimensions:\n")
print(dim(A_example))

cat("\nExample concordance values:\n")
print(a_quad(A_example))


# ============================================================
# STEP 2: Monte Carlo approximation to tau*
# ============================================================

tau_star <- function(X, Y, M = 50000, verbose = FALSE) {
  
  n <- nrow(X)
  
  if (n < 4) {
    stop("At least four observations are required.")
  }
  
  if (verbose) {
    cat("\n--- Estimating tau* ---\n")
    cat("Sample size:", n, "\n")
    cat("Number of Monte Carlo quadruples:", M, "\n")
  }
  
  ids <- matrix(
    sample.int(n, 4 * M, replace = TRUE),
    nrow = 4
  )
  
  bad <- apply(ids, 2, function(z) {
    length(unique(z)) < 4
  })
  
  while (any(bad)) {
    
    ids[, bad] <- matrix(
      sample.int(n, 4 * sum(bad), replace = TRUE),
      nrow = 4
    )
    
    bad <- apply(ids, 2, function(z) {
      length(unique(z)) < 4
    })
  }
  
  if (verbose) {
    cat("Valid quadruples generated:", ncol(ids), "\n")
    cat("First five quadruples:\n")
    print(ids[, seq_len(min(5, ncol(ids))), drop = FALSE])
  }
  
  Xarr <- array(
    NA_real_,
    dim = c(4, ncol(X), M)
  )
  
  Yarr <- array(
    NA_real_,
    dim = c(4, ncol(Y), M)
  )
  
  for (k in 1:4) {
    
    Xarr[k, , ] <- t(X[ids[k, ], , drop = FALSE])
    Yarr[k, , ] <- t(Y[ids[k, ], , drop = FALSE])
  }
  
  concordance_X <- a_quad(Xarr)
  concordance_Y <- a_quad(Yarr)
  
  products <- concordance_X * concordance_Y
  
  tau_hat <- mean(products)
  
  if (verbose) {
    
    cat("\nConcordance summary for X:\n")
    print(summary(concordance_X))
    
    cat("\nConcordance summary for Y:\n")
    print(summary(concordance_Y))
    
    cat("\nProduct summary:\n")
    print(summary(products))
    
    cat("\nEstimated tau*:\n")
    print(tau_hat)
  }
  
  tau_hat
}


# ============================================================
# STEP 3: Estimate the projection variance
# ============================================================

estimate_variance <- function(
    X, Y, Mvar = 20000, verbose = FALSE) {
  
  n <- nrow(X)
  
  if (verbose) {
    cat("\n--- Estimating projection variance ---\n")
    cat("Sample size:", n, "\n")
    cat("Requested Monte Carlo draws:", Mvar, "\n")
  }
  
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
    
    xarr <- array(
      NA_real_,
      dim = c(4, ncol(X), 1)
    )
    
    yarr <- array(
      NA_real_,
      dim = c(4, ncol(Y), 1)
    )
    
    for (j in 1:4) {
      xarr[j, , 1] <- X[ids[j], ]
      yarr[j, , 1] <- Y[ids[j], ]
    }
    
    h2[k] <- a_quad(xarr) * a_quad(yarr)
    
    # Display occasional intermediate values
    if (verbose && k %% 5000 == 0) {
      cat(
        "Variance simulation:",
        k, "of", Mvar,
        "| Current valid draws:",
        sum(is.finite(h2[seq_len(k)])),
        "\n"
      )
    }
  }
  
  h2 <- h2[is.finite(h2)]
  
  if (verbose) {
    cat("\nValid variance draws:", length(h2), "\n")
    cat("Invalid draws discarded:", Mvar - length(h2), "\n")
    cat("\nProjection-kernel summary:\n")
    print(summary(h2))
  }
  
  if (length(h2) < 2) {
    return(NA_real_)
  }
  
  sigma2_hat <- 4 * var(h2)
  
  if (verbose) {
    cat("\nEstimated variance (sigma squared):\n")
    print(sigma2_hat)
  }
  
  sigma2_hat
}


# ============================================================
# STEP 4: Independence test
# ============================================================

independence_test <- function(
    n, B = 1000, M = 50000, Mvar = 20000) {
  
  cat("\n\n")
  cat("##################################################\n")
  cat("RUNNING INDEPENDENCE TEST FOR n =", n, "\n")
  cat("##################################################\n")
  
  # ----------------------------------------------------------
  # 4.1 Generate correlated bivariate normal observations
  # ----------------------------------------------------------
  
  cat("\nSTEP 4.1: Generating the dataset\n")
  
  Sigma <- matrix(
    c(1, rho, rho, 1),
    nrow = 2
  )
  
  cat("Correlation matrix:\n")
  print(Sigma)
  
  dat <- MASS::mvrnorm(
    n = n,
    mu = c(0, 0),
    Sigma = Sigma
  )
  
  X <- matrix(dat[, 1], ncol = 1)
  Y <- matrix(dat[, 2], ncol = 1)
  
  cat("\nFirst six generated observations:\n")
  
  print(
    head(
      data.frame(X = X[, 1], Y = Y[, 1]),
      6
    )
  )
  
  cat("\nDimensions of X and Y:\n")
  print(dim(X))
  print(dim(Y))
  
  cat("\nSummary of X:\n")
  print(summary(X))
  
  cat("\nSummary of Y:\n")
  print(summary(Y))
  
  cat("\nObserved Pearson correlation:\n")
  print(cor(X[, 1], Y[, 1]))
  
  # ----------------------------------------------------------
  # 4.2 Estimate tau*
  # ----------------------------------------------------------
  
  cat("\nSTEP 4.2: Calculating observed tau*\n")
  
  tau_obs <- tau_star(
    X, Y,
    M = M,
    verbose = show_details
  )
  
  cat("\nObserved tau* =", tau_obs, "\n")
  
  # ----------------------------------------------------------
  # 4.3 Estimate variance
  # ----------------------------------------------------------
  
  cat("\nSTEP 4.3: Estimating variance\n")
  
  sigma2_hat <- estimate_variance(
    X, Y,
    Mvar = Mvar,
    verbose = show_details
  )
  
  cat("\nEstimated sigma squared =", sigma2_hat, "\n")
  
  # ----------------------------------------------------------
  # 4.4 Asymptotic test
  # ----------------------------------------------------------
  
  cat("\nSTEP 4.4: Asymptotic test\n")
  
  if (is.finite(sigma2_hat) && sigma2_hat > 0) {
    
    z_stat <- sqrt(n) * tau_obs / sqrt(sigma2_hat)
    
    p_asym <- 2 * pnorm(-abs(z_stat))
    
  } else {
    
    z_stat <- NA_real_
    p_asym <- NA_real_
  }
  
  asym_decision <- if (is.na(p_asym)) {
    NA_character_
  } else if (p_asym <= alpha) {
    "Reject H0"
  } else {
    "Do not reject H0"
  }
  
  cat("Z statistic:\n")
  print(z_stat)
  
  cat("Asymptotic p-value:\n")
  print(p_asym)
  
  cat("Asymptotic decision:\n")
  print(asym_decision)
  
  # ----------------------------------------------------------
  # 4.5 Permutation test
  # ----------------------------------------------------------
  
  cat("\nSTEP 4.5: Permutation test\n")
  cat("Total permutations:", B, "\n")
  
  tau_perm <- numeric(B)
  
  for (b in seq_len(B)) {
    
    Y_perm <- Y[sample.int(n), , drop = FALSE]
    
    tau_perm[b] <- tau_star(
      X, Y_perm,
      M = M,
      verbose = FALSE
    )
    
    # Display progress every 50 permutations
    if (b %% 50 == 0 || b == B) {
      
      cat(
        "Permutation:", b, "/", B,
        "| Latest tau*:", round(tau_perm[b], 5),
        "| Current permutation p-value:",
        round(
          (1 + sum(abs(tau_perm[seq_len(b)]) >= abs(tau_obs))) /
            (b + 1),
          5
        ),
        "\n"
      )
    }
  }
  
  p_perm <- (
    1 + sum(abs(tau_perm) >= abs(tau_obs))
  ) / (B + 1)
  
  perm_decision <- if (p_perm <= alpha) {
    "Reject H0"
  } else {
    "Do not reject H0"
  }
  
  cat("\nSummary of permutation tau* values:\n")
  print(summary(tau_perm))
  
  cat("\nFirst ten permutation statistics:\n")
  print(head(tau_perm, 10))
  
  cat("\nPermutation p-value:\n")
  print(p_perm)
  
  cat("\nPermutation decision:\n")
  print(perm_decision)
  
  # ----------------------------------------------------------
  # 4.6 Return results
  # ----------------------------------------------------------
  
  result <- data.frame(
    Sample_Size = n,
    Tau_Star = tau_obs,
    Sigma_Squared = sigma2_hat,
    Z_Statistic = z_stat,
    Asymptotic_p = p_asym,
    Permutation_p = p_perm,
    Asymptotic_Decision = asym_decision,
    Permutation_Decision = perm_decision
  )
  
  cat("\nFINAL RESULT FOR n =", n, "\n")
  print(result, row.names = FALSE)
  
  result
}


# ============================================================
# STEP 5: Run the test for all sample sizes
# ============================================================

cat("\n\n")
cat("====================================================\n")
cat("STARTING COMPLETE SIMULATION\n")
cat("Sample sizes:", sample_sizes, "\n")
cat("====================================================\n")

results_list <- vector("list", length(sample_sizes))

for (i in seq_along(sample_sizes)) {
  
  cat("\n\n")
  cat("************ SAMPLE SIZE INDEX", i, "************\n")
  
  results_list[[i]] <- independence_test(
    n = sample_sizes[i],
    B = B,
    M = M,
    Mvar = Mvar
  )
  
  cat("\nResults obtained so far:\n")
  
  interim_results <- do.call(rbind, results_list[seq_len(i)])
  
  interim_results[, 2:6] <- round(
    interim_results[, 2:6],
    4
  )
  
  print(interim_results, row.names = FALSE)
}


# ============================================================
# STEP 6: Final results table
# ============================================================

results <- do.call(rbind, results_list)

results[, 2:6] <- round(
  results[, 2:6],
  4
)

cat("\n\n")
cat("====================================================\n")
cat("FINAL NUMERICAL RESULTS\n")
cat("====================================================\n")

print(results, row.names = FALSE)


# ============================================================
# STEP 7: Save final results
# ============================================================

write.csv(
  results,
  file = "independence_test_results.csv",
  row.names = FALSE
)

cat("\nFinal results saved to:\n")
cat(file.path(getwd(), "independence_test_results.csv"), "\n")

cat("\nSIMULATION COMPLETED SUCCESSFULLY\n")
