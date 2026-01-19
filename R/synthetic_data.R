# Helper function to generate correlated FSC-A and FSC-H for singlets
generate_singlets <- function(n, fsc_mean, fsc_sd) {
  fsc_a <- rlnorm(n, log(fsc_mean), fsc_sd)
  # FSC-H should be proportional to FSC-A for singlets with some noise
  fsc_h <- fsc_a * rnorm(n, mean = 0.85, sd = 0.05)
  list(fsc_a = fsc_a, fsc_h = fsc_h)
}

# Helper function to generate doublets (higher FSC-A for same FSC-H)
generate_doublets <- function(n, fsc_mean, fsc_sd) {
  fsc_h <- rlnorm(n, log(fsc_mean * 0.7), fsc_sd)
  fsc_a <- fsc_h * rnorm(n, mean = 1.5, sd = 0.1)  # Higher FSC-A
  list(fsc_a = fsc_a, fsc_h = fsc_h)
}

# Helper function to generate debris (low FSC/SSC)
generate_debris <- function(n) {
  fsc_a <- rlnorm(n, log(5000), 0.5)
  ssc_a <- rlnorm(n, log(3000), 0.5)
  fsc_h <- fsc_a * rnorm(n, mean = 0.7, sd = 0.2)
  list(fsc_a = fsc_a, ssc_a = ssc_a, fsc_h = fsc_h)
}

# Helper function to generate live cells (low FL1-A)
generate_live_cells <- function(n) {
  rlnorm(n, 0, 0.5)
}

# Helper function to generate dead cells (high FL1-A)
generate_dead_cells <- function(n) {
  rlnorm(n, log(50000), 0.3)
}

# Helper function to generate marker populations
generate_marker <- function(n, pct_positive) {
  n_pos <- round(n * pct_positive)
  n_neg <- n - n_pos

  # Negative population (low fluorescence)
  neg <- rlnorm(n_neg, log(300), 0.3)
  # Positive population (high fluorescence)
  pos <- rlnorm(n_pos, log(15000), 0.4)

  sample(c(neg, pos))  # Shuffle
}