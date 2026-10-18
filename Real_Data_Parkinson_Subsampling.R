############################################################
# REAL DATA APPLICATION:
# PARKINSON'S TELEMONITORING DATASET
#
# Objective:
# Examine dependence between multivariate voice biomarkers
# and motor_UPDRS using:
#   1. Custom multivariate sign-product statistic
#   2. Permutation test
#   3. Distance correlation
#   4. Maximal Information Coefficient (MIC)
#   5. Subsampling-based confidence interval
#
# IMPORTANT:
# The custom sign-product statistic below reproduces the
# construction in the existing manuscript code. It is not
# automatically equivalent to canonical Bergsma-Dassios tau*.
############################################################


############################################################
# 1. INSTALL AND LOAD REQUIRED PACKAGES
############################################################

required_packages <- c("energy", "minerva", "ggplot2",
                       "dplyr", "tidyr", "readr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

library(energy)
library(minerva)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)


############################################################
# 2. SET REPRODUCIBILITY AND COMPUTATIONAL SETTINGS
############################################################

set.seed(123)

# Number of observations used in the computational analysis.
# The complete dataset is retained for descriptive analysis.
# Increase this value after checking computational resources.
analysis_n <- 120

# Number of random quadruples used to approximate the
# custom sign-product statistic.
M_quadruples <- 10000

# Number of permutations for hypothesis testing.
B_permutation <- 999

# Number of subsamples for subsampling inference.
M_subsamples <- 500

# Subsample-size exponent:
# b_n = floor(n^gamma)
gamma <- 0.7

alpha <- 0.05


############################################################
# 3. EXTRACT AND READ THE UPLOADED DATASET
############################################################

# 1. Define paths
zip_file <- "C:/Users/hp/Desktop/parkinsons+telemonitoring.zip"
output_dir <- "C:/Users/hp/Documents/Parkinson-Output-Values"

# 2. Extract all contents
unzip(zip_file, exdir = output_dir)

print(zip_file)

data_member <- "parkinsons_updrs.data"
# Extract the data file into a temporary directory.
extract_dir <- tempfile("parkinsons_data_")
dir.create(extract_dir)

unzip(zip_file, files = data_member, exdir = extract_dir)

data_path <- file.path(extract_dir, data_member)

# Preserve the original UCI column names, including punctuation.
dat <- read.csv(
  data_path,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


############################################################
# 4. INITIAL DATA INSPECTION
############################################################

cat("\n---------------- DATA DIMENSIONS ----------------\n")
print(dim(dat))

cat("\n---------------- COLUMN NAMES -------------------\n")
print(names(dat))

cat("\n---------------- FIRST SIX OBSERVATIONS ---------\n")
print(head(dat))

cat("\n---------------- DATA STRUCTURE -----------------\n")
str(dat)

cat("\n---------------- MISSING VALUES -----------------\n")
print(colSums(is.na(dat)))

cat("\n---------------- DUPLICATE ROWS -----------------\n")
print(sum(duplicated(dat)))


############################################################
# 5. VALIDATE REQUIRED VARIABLES
############################################################

required_variables <- c(
  "subject#",
  "age",
  "sex",
  "test_time",
  "motor_UPDRS",
  "total_UPDRS",
  "Jitter(%)",
  "Shimmer",
  "NHR",
  "HNR"
)

missing_variables <- setdiff(required_variables, names(dat))

if (length(missing_variables) > 0) {
  stop(
    paste(
      "The following required variables are missing:",
      paste(missing_variables, collapse = ", ")
    )
  )
}


############################################################
# 6. DEFINE THE ANALYTIC VARIABLES
############################################################

# Voice biomarkers:
# Jitter(%), Shimmer, NHR, and HNR.
#
# Clinical outcome:
# motor_UPDRS.

voice_variables <- c(
  "Jitter(%)",
  "Shimmer",
  "NHR",
  "HNR"
)

outcome_variable <- "motor_UPDRS"

analysis_variables <- c(
  voice_variables,
  outcome_variable
)

# Select the required variables.
dat_analysis <- dat[, analysis_variables, drop = FALSE]

# Convert to numeric if necessary.
dat_analysis[] <- lapply(
  dat_analysis,
  function(x) as.numeric(as.character(x))
)

# Remove rows containing missing or non-finite values
# in the selected variables.
valid_rows <- complete.cases(dat_analysis)

for (v in names(dat_analysis)) {
  valid_rows <- valid_rows & is.finite(dat_analysis[[v]])
}

dat_analysis <- dat_analysis[valid_rows, , drop = FALSE]

cat("\nNumber of rows in original dataset:", nrow(dat), "\n")
cat("Number of complete analytic rows:", nrow(dat_analysis), "\n")

if (nrow(dat_analysis) < 30) {
  stop("Too few complete observations for this analysis.")
}


############################################################
# 7. DESCRIPTIVE STATISTICS
############################################################

descriptive_statistics <- data.frame(
  Variable = names(dat_analysis),
  N = sapply(dat_analysis, length),
  Mean = sapply(dat_analysis, mean),
  SD = sapply(dat_analysis, sd),
  Minimum = sapply(dat_analysis, min),
  Median = sapply(dat_analysis, median),
  Maximum = sapply(dat_analysis, max),
  row.names = NULL
)

descriptive_statistics[, -1] <-
  round(descriptive_statistics[, -1], 5)

cat("\n---------------- DESCRIPTIVE STATISTICS ---------\n")
print(descriptive_statistics)

write.csv(
  descriptive_statistics,
  "parkinsons_descriptive_statistics.csv",
  row.names = FALSE
)


############################################################
# 8. SUBJECT-LEVEL INFORMATION
############################################################

# Repeated voice recordings are available for each subject.
# Therefore, observations should not automatically be
# interpreted as independent patients.
# Summarise recordings by subject using base R
subject_summary <- do.call(rbind,lapply(split(dat, dat[["subject#"]]), function(d) {
      x <- d$motor_UPDRS
        data.frame(
          subject = d[["subject#"]][1],
          Number_of_recordings = nrow(d),
          Mean_motor_UPDRS = if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE),
          SD_motor_UPDRS = if (sum(!is.na(x)) < 2) NA_real_ else sd(x, na.rm = TRUE)
          )
      })
    )
# Display the results
print(subject_summary)

cat("\n---------------- SUBJECT SUMMARY ----------------\n")
print(summary(subject_summary$Number_of_recordings))

write.csv(
  subject_summary,
  "parkinsons_subject_summary.csv",
  row.names = FALSE
)


############################################################
# 9. EXPLORATORY PLOTS
############################################################

############################################################
# 9.1 HISTOGRAMS OF SELECTED VARIABLES
############################################################

# Convert selected analysis variables to numeric
dat_hist <- as.data.frame(dat_analysis)

for (j in seq_along(dat_hist)) {
  if (is.factor(dat_hist[[j]])) {
    dat_hist[[j]] <- as.numeric(as.character(dat_hist[[j]]))
  } else if (is.character(dat_hist[[j]])) {
    dat_hist[[j]] <- as.numeric(dat_hist[[j]])
  }
}

# Keep numeric columns only
numeric_columns <- names(dat_hist)[
  sapply(dat_hist, is.numeric)
]

if (length(numeric_columns) == 0) {
  stop("Error: dat_analysis contains no numeric columns.")
}

dat_hist <- dat_hist[, numeric_columns, drop = FALSE]

# Reshape data using base R (no pivot_longer or pipes)
long_data <- do.call(
  rbind,
  lapply(names(dat_hist), function(v) {
    data.frame(
      Variable = v,
      Value = dat_hist[[v]]
    )
  })
)

# Remove missing and non-finite values
long_data <- long_data[
  is.finite(long_data$Value),
  ,
  drop = FALSE
]

# Plot histograms
p_hist <- ggplot(
  long_data,
  aes(x = Value)
) +
  geom_histogram(
    bins = 35,
    color = "black",
    fill = "lightblue"
  ) +
  facet_wrap(
    ~ Variable,
    scales = "free",
    ncol = 2
  ) +
  theme_bw() +
  labs(
    title = "Distributions of Voice Biomarkers and Motor UPDRS",
    x = "Observed value",
    y = "Frequency"
  )

print(p_hist)

ggsave(
  filename = "parkinsons_variable_histograms.png",
  plot = p_hist,
  width = 9,
  height = 7,
  dpi = 300
)


############################################################
# 9.2 SCATTERPLOTS OF VOICE BIOMARKERS AGAINST MOTOR UPDRS
############################################################

# Check required variables
if (!exists("voice_variables")) {
  stop("Error: voice_variables has not been defined.")
}

if (!("motor_UPDRS" %in% names(dat_analysis))) {
  stop("Error: motor_UPDRS is not a column in dat_analysis.")
}

missing_variables <- setdiff(
  voice_variables,
  names(dat_analysis)
)

if (length(missing_variables) > 0) {
  stop(
    paste(
      "These voice_variables are missing from dat_analysis:",
      paste(missing_variables, collapse = ", ")
    )
  )
}

# Select the required columns
scatter_data <- dat_analysis[
  ,
  c(voice_variables, "motor_UPDRS"),
  drop = FALSE
]

# Convert selected columns to numeric if necessary
for (v in names(scatter_data)) {
  if (is.factor(scatter_data[[v]]) ||
      is.character(scatter_data[[v]])) {
    scatter_data[[v]] <- as.numeric(
      as.character(scatter_data[[v]])
    )
  }
}

# Reshape the biomarker columns using base R
scatter_data <- do.call(
  rbind,
  lapply(voice_variables, function(v) {
    data.frame(
      Biomarker = v,
      Biomarker_Value = scatter_data[[v]],
      motor_UPDRS = scatter_data[["motor_UPDRS"]]
    )
  })
)

# Remove missing and non-finite observations
scatter_data <- scatter_data[
  is.finite(scatter_data$Biomarker_Value) &
    is.finite(scatter_data$motor_UPDRS),
  ,
  drop = FALSE
]

# Create scatterplots with LOESS curves
p_scatter <- ggplot(
  scatter_data,
  aes(
    x = Biomarker_Value,
    y = motor_UPDRS
  )
) +
  geom_point(alpha = 0.35) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = TRUE
  ) +
  facet_wrap(
    ~ Biomarker,
    scales = "free_x",
    ncol = 2
  ) +
  theme_bw() +
  labs(
    title = "Voice Biomarkers versus Motor UPDRS",
    x = "Voice biomarker",
    y = "Motor UPDRS"
  )

print(p_scatter)

ggsave(
  filename = "parkinsons_voice_vs_motor_UPDRS.png",
  plot = p_scatter,
  width = 10,
  height = 7,
  dpi = 300
)

############################################################
# 10. STANDARDIZATION
############################################################

# Standardization improves numerical comparability.
# Rank/sign-based comparisons are invariant to strictly
# increasing affine transformations of individual variables,
# apart from possible effects of ties and numerical precision.

X_all <- as.matrix(dat_analysis[, voice_variables])

Y_all <- as.numeric(dat_analysis[[outcome_variable]])

X_all <- scale(X_all)
Y_all <- as.numeric(scale(Y_all))

if (any(!is.finite(X_all)) || any(!is.finite(Y_all))) {
  stop("Non-finite values encountered after standardization.")
}


############################################################
# 11. DRAW A COMMON ANALYSIS SUBSAMPLE
############################################################

# For computational comparability, the proposed statistic,
# distance correlation, and MIC are evaluated on the same
# randomly selected observations.

set.seed(123)

n_available <- nrow(X_all)
n_analysis <- min(analysis_n, n_available)

analysis_index <- sample(
  seq_len(n_available),
  size = n_analysis,
  replace = FALSE
)

X <- X_all[analysis_index, , drop = FALSE]
Y <- Y_all[analysis_index]

n <- nrow(X)

cat("\n---------------- ANALYSIS SAMPLE ----------------\n")
cat("Analysis sample size:", n, "\n")
cat("Number of voice biomarkers:", ncol(X), "\n")

b_n <- floor(n^gamma)

if (b_n < 4 || b_n >= n) {
  stop("Invalid subsample size. Check analysis_n and gamma.")
}

cat("Subsample size b_n:", b_n, "\n")


############################################################
# 12. DEFINE THE CUSTOM SIGN-PRODUCT STATISTIC
############################################################

# This function follows the construction in the existing
# manuscript code:
#
# sign.vec(x1,x2) = product_j sign(x1j - x2j)
#
# For each selected quadruple (i,j,l,r), the kernel is
#
# sign.vec(X_i,X_j) * sign.vec(X_l,X_r) *
# sign(Y_i-Y_j) * sign(Y_l-Y_r).
#
# WARNING:
# This is a custom sign-product statistic. Its equivalence
# to canonical multivariate Bergsma-Dassios tau* must not
# be assumed without a mathematical proof.

sign_vec <- function(x1, x2) {
  prod(sign(x1 - x2))
}


# Generate random quadruples once, and reuse them when
# evaluating the statistic on the same sample.
#
# Reusing the quadruples during permutation reduces
# additional Monte Carlo noise from repeatedly drawing
# different quadruples.

make_quadruples <- function(n, M) {
  
  if (n < 4) {
    stop("At least four observations are required.")
  }
  
  matrix(
    sample.int(n, size = 4 * M, replace = TRUE),
    ncol = 4
  )
}


# Evaluate the custom statistic on a given sample and a
# supplied matrix of quadruple indices.
#
# Repeated indices in a quadruple are discarded because
# the statistic is intended to use four distinct observations.

custom_tau_from_quadruples <- function(X, Y, quadruples) {
  
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  
  n <- nrow(X)
  
  if (length(Y) != n) {
    stop("X and Y must contain the same number of observations.")
  }
  
  if (ncol(quadruples) != 4) {
    stop("The quadruple index matrix must have four columns.")
  }
  
  valid <- apply(
    quadruples,
    1,
    function(z) length(unique(z)) == 4
  )
  
  quadruples <- quadruples[valid, , drop = FALSE]
  
  if (nrow(quadruples) == 0) {
    stop("No valid quadruples were generated.")
  }
  
  kernel_values <- numeric(nrow(quadruples))
  
  for (k in seq_len(nrow(quadruples))) {
    
    i <- quadruples[k, 1]
    j <- quadruples[k, 2]
    l <- quadruples[k, 3]
    r <- quadruples[k, 4]
    
    sx1 <- sign_vec(X[i, ], X[j, ])
    sx2 <- sign_vec(X[l, ], X[r, ])
    
    sy1 <- sign(Y[i] - Y[j])
    sy2 <- sign(Y[l] - Y[r])
    
    kernel_values[k] <- sx1 * sx2 * sy1 * sy2
  }
  
  mean(kernel_values)
}


# Draw quadruples and compute the statistic.
set.seed(123)

quadruples <- make_quadruples(
  n = n,
  M = M_quadruples
)

tau_observed <- custom_tau_from_quadruples(
  X = X,
  Y = Y,
  quadruples = quadruples
)

cat("\n---------------- CUSTOM SIGN-PRODUCT STATISTIC ----------------\n")
cat("Observed statistic:", tau_observed, "\n")
cat("Number of sampled quadruples:", nrow(quadruples), "\n")


############################################################
# 13. PERMUTATION TEST FOR THE CUSTOM STATISTIC
############################################################

# Null hypothesis:
# H0: X and Y are independent.
#
# Under the null, the Y values are permuted relative to X.
#
# The validity of this permutation scheme assumes exchangeable
# observations under the null. Because this dataset contains
# repeated recordings within subjects, unrestricted row-wise
# permutation may not be valid for population-level inference.
# A subject-aware procedure may be necessary.

set.seed(456)

tau_permuted <- numeric(B_permutation)

for (b in seq_len(B_permutation)) {
  
  perm_index <- sample.int(n)
  
  tau_permuted[b] <- custom_tau_from_quadruples(
    X = X,
    Y = Y[perm_index],
    quadruples = quadruples
  )
  
  if (b %% 100 == 0) {
    cat("Completed permutation:", b, "\n")
  }
}


# Two-sided permutation p-value.
p_tau <- (
  1 + sum(abs(tau_permuted) >= abs(tau_observed))
) / (B_permutation + 1)

cat("\n---------------- PERMUTATION TEST ----------------\n")
cat("Observed statistic:", tau_observed, "\n")
cat("Number of permutations:", B_permutation, "\n")
cat("Two-sided permutation p-value:", p_tau, "\n")


# Plot the permutation distribution.
perm_data <- data.frame(Statistic = tau_permuted)

p_perm <- ggplot(perm_data, aes(x = Statistic)) +
  geom_histogram(
    bins = 35,
    color = "black"
  ) +
  geom_vline(
    xintercept = tau_observed,
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_vline(
    xintercept = -abs(tau_observed),
    linetype = "dotted"
  ) +
  geom_vline(
    xintercept = abs(tau_observed),
    linetype = "dotted"
  ) +
  theme_bw() +
  labs(
    title = "Permutation Distribution of the Custom Statistic",
    subtitle = paste("Permutation p-value =", signif(p_tau, 4)),
    x = "Permuted statistic",
    y = "Frequency"
  )

ggsave(
  "parkinsons_custom_statistic_permutation.png",
  p_perm,
  width = 8,
  height = 5,
  dpi = 300
)


############################################################
# 14. DISTANCE CORRELATION
############################################################

# energy::dcor computes distance correlation.
# Here, X is a multivariate matrix and Y is a univariate vector.

dcor_observed <- energy::dcor(X, Y)

cat("\n---------------- DISTANCE CORRELATION ----------------\n")
cat("Observed distance correlation:", dcor_observed, "\n")


# Permutation test for distance correlation.
set.seed(789)

dcor_permuted <- numeric(B_permutation)

for (b in seq_len(B_permutation)) {
  
  perm_index <- sample.int(n)
  
  dcor_permuted[b] <- energy::dcor(
    X,
    Y[perm_index]
  )
  
  if (b %% 100 == 0) {
    cat("Distance-correlation permutation:", b, "\n")
  }
}

p_dcor <- (
  1 + sum(dcor_permuted >= dcor_observed)
) / (B_permutation + 1)

cat("Distance correlation:", dcor_observed, "\n")
cat("Distance correlation permutation p-value:", p_dcor, "\n")


############################################################
# 15. MAXIMAL INFORMATION COEFFICIENT (MIC)
############################################################

# MIC is fundamentally a bivariate association measure.
# Calculate MIC for each voice biomarker against motor_UPDRS.
#
# These individual MIC values are not a single multivariate
# test of independence between the full vector X and Y.

mic_results <- data.frame(
  Biomarker = voice_variables,
  MIC = NA_real_,
  stringsAsFactors = FALSE
)

for (j in seq_along(voice_variables)) {
  
  mic_fit <- minerva::mine(
    cbind(
      X[, j],
      Y
    ),
    use = "pairwise.complete.obs"
  )
  
  mic_results$MIC[j] <- as.numeric(mic_fit$MIC[1, 2])
}

cat("\n---------------- MIC RESULTS ----------------\n")
print(mic_results)

write.csv(
  mic_results,
  "parkinsons_MIC_results.csv",
  row.names = FALSE
)


############################################################
# 16. PERMUTATION TESTS FOR INDIVIDUAL MIC VALUES
############################################################

# For each biomarker, permute Y and compare the observed MIC
# with the permutation distribution.
#
# These are four separate bivariate tests. Apply a multiple-
# testing correction when interpreting them together.

set.seed(987)

mic_test_results <- data.frame(
  Biomarker = voice_variables,
  Observed_MIC = NA_real_,
  P_value = NA_real_,
  stringsAsFactors = FALSE
)

for (j in seq_along(voice_variables)) {
  
  xj <- X[, j]
  
  observed_fit <- minerva::mine(
    cbind(xj, Y),
    use = "pairwise.complete.obs"
  )
  
  observed_mic <- as.numeric(observed_fit$MIC[1, 2])
  
  mic_perm <- numeric(B_permutation)
  
  for (b in seq_len(B_permutation)) {
    
    perm_index <- sample.int(n)
    
    perm_fit <- minerva::mine(
      cbind(xj, Y[perm_index]),
      use = "pairwise.complete.obs"
    )
    
    mic_perm[b] <- as.numeric(perm_fit$MIC[1, 2])
  }
  
  p_mic <- (
    1 + sum(mic_perm >= observed_mic)
  ) / (B_permutation + 1)
  
  mic_test_results$Observed_MIC[j] <- observed_mic
  mic_test_results$P_value[j] <- p_mic
  
  cat(
    "\nMIC test for",
    voice_variables[j],
    ": MIC =",
    observed_mic,
    ", p-value =",
    p_mic,
    "\n"
  )
}

mic_test_results$P_value_Holm <- p.adjust(
  mic_test_results$P_value,
  method = "holm"
)

print(mic_test_results)

write.csv(
  mic_test_results,
  "parkinsons_MIC_permutation_tests.csv",
  row.names = FALSE
)


############################################################
# 17. SUBSAMPLING-BASED INFERENCE FOR THE CUSTOM STATISTIC
############################################################

# The following implements the centered subsampling approach
# in the existing manuscript code.
#
# It is an exploratory interval for the custom statistic.
# It should not be presented as validated inference until
# the statistic's asymptotic rate and subsampling theory
# have been justified.

set.seed(2468)

b_n <- floor(n^gamma)

subsample_statistics <- numeric(M_subsamples)

for (m in seq_len(M_subsamples)) {
  
  idx <- sample.int(
    n,
    size = b_n,
    replace = FALSE
  )
  
  # Each subsample gets its own random quadruples.
  # Increase M_quadruples if more accurate evaluation is needed.
  q_sub <- make_quadruples(
    n = b_n,
    M = M_quadruples
  )
  
  subsample_statistics[m] <- custom_tau_from_quadruples(
    X = X[idx, , drop = FALSE],
    Y = Y[idx],
    quadruples = q_sub
  )
  
  if (m %% 50 == 0) {
    cat("Completed subsample:", m, "\n")
  }
}


# Center and scale the subsample estimates.
subsampling_root_b <- sqrt(b_n) * (
  subsample_statistics - tau_observed
)

# Quantiles of the centered subsampling distribution.
q_lower <- unname(
  quantile(subsampling_root_b, probs = alpha / 2)
)

q_upper <- unname(
  quantile(subsampling_root_b, probs = 1 - alpha / 2)
)


# Candidate percentile-type interval using the root-n scaling
# formula in the manuscript.
#
# This interval is exploratory until the required asymptotic
# conditions have been established.

ci_lower <- tau_observed - q_upper / sqrt(n)
ci_upper <- tau_observed - q_lower / sqrt(n)

cat("\n---------------- SUBSAMPLING RESULTS ----------------\n")
cat("Full-sample custom statistic:", tau_observed, "\n")
cat("Subsample size:", b_n, "\n")
cat("Number of subsamples:", M_subsamples, "\n")
cat("Nominal confidence level:", 1 - alpha, "\n")
cat("Candidate lower confidence limit:", ci_lower, "\n")
cat("Candidate upper confidence limit:", ci_upper, "\n")
cat("Candidate interval width:", ci_upper - ci_lower, "\n")


############################################################
# 18. PLOT THE SUBSAMPLING DISTRIBUTION
############################################################

subsampling_data <- data.frame(
  Statistic = subsample_statistics
)

p_subsampling <- ggplot(
  subsampling_data,
  aes(x = Statistic)
) +
  geom_histogram(
    bins = 35,
    color = "black"
  ) +
  geom_vline(
    xintercept = tau_observed,
    linetype = "dashed",
    linewidth = 1
  ) +
  theme_bw() +
  labs(
    title = "Subsampling Distribution of the Custom Statistic",
    subtitle = paste(
      "Subsample size =", b_n,
      "; number of subsamples =", M_subsamples
    ),
    x = "Subsample statistic",
    y = "Frequency"
  )

ggsave(
  "parkinsons_subsampling_distribution.png",
  p_subsampling,
  width = 8,
  height = 5,
  dpi = 300
)


############################################################
# 19. SUMMARIZE ALL DEPENDENCE RESULTS
############################################################

comparison_results <- data.frame(
  Method = c(
    "Custom sign-product statistic",
    "Distance correlation"
  ),
  Statistic = c(
    tau_observed,
    dcor_observed
  ),
  P_value = c(
    p_tau,
    p_dcor
  )
)

cat("\n---------------- METHOD COMPARISON ----------------\n")
print(comparison_results)

write.csv(
  comparison_results,
  "parkinsons_dependence_comparison.csv",
  row.names = FALSE
)


############################################################
# 20. SAVE ALL MAIN NUMERICAL RESULTS
############################################################

main_results <- data.frame(
  Method = c(
    "Custom sign-product statistic",
    "Distance correlation"
  ),
  Statistic = c(
    tau_observed,
    dcor_observed
  ),
  P_value = c(
    p_tau,
    p_dcor
  ),
  Sample_Size = n,
  Permutations = B_permutation,
  stringsAsFactors = FALSE
)

write.csv(
  main_results,
  "parkinsons_main_results.csv",
  row.names = FALSE
)


############################################################
# 21. SAVE THE ANALYSIS SETTINGS
############################################################

analysis_settings <- data.frame(
  Setting = c(
    "Number of complete analytic rows",
    "Analysis sample size",
    "Number of voice biomarkers",
    "Number of sampled quadruples",
    "Number of permutations",
    "Number of subsamples",
    "Subsample-size exponent",
    "Subsample size",
    "Nominal significance level",
    "Random seed"
  ),
  Value = c(
    nrow(dat_analysis),
    n,
    ncol(X),
    M_quadruples,
    B_permutation,
    M_subsamples,
    gamma,
    b_n,
    alpha,
    123
  )
)

write.csv(
  analysis_settings,
  "parkinsons_analysis_settings.csv",
  row.names = FALSE
)


############################################################
# 22. FINAL CONSOLE SUMMARY
############################################################

cat("\n====================================================\n")
cat("PARKINSON'S REAL DATA ANALYSIS COMPLETED\n")
cat("====================================================\n")

cat("\nCustom sign-product statistic:", tau_observed)
cat("\nCustom statistic permutation p-value:", p_tau)

cat("\n\nDistance correlation:", dcor_observed)
cat("\nDistance correlation permutation p-value:", p_dcor)

cat("\n\nCandidate subsampling confidence interval: [",
    ci_lower, ", ", ci_upper, "]", sep = "")

cat("\n\nFiles written to the current working directory:\n")

output_files <- c(
  "parkinsons_descriptive_statistics.csv",
  "parkinsons_subject_summary.csv",
  "parkinsons_MIC_results.csv",
  "parkinsons_MIC_permutation_tests.csv",
  "parkinsons_dependence_comparison.csv",
  "parkinsons_main_results.csv",
  "parkinsons_analysis_settings.csv",
  "parkinsons_variable_histograms.png",
  "parkinsons_voice_vs_motor_UPDRS.png",
  "parkinsons_custom_statistic_permutation.png",
  "parkinsons_subsampling_distribution.png"
)

print(output_files)

cat("\n====================================================\n")
