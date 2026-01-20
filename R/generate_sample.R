generate_sample <- function(
  n_events = 20000,
  debris_pct = runif(1, 0.01, 0.15),
  live_pct = runif(1, 0.75, 0.99),
  size_pct = runif(1, 0.75, 0.99),
  singlet_pct = runif(1, 0.85, 0.99),
  fl2_positive_pct = runif(1, 0.05, 0.95),
  fl3_positive_pct = runif(1, 0.05, 0.95),
  fsc_mean = runif(1, 1e5, 2e5),
  ssc_mean = runif(1, 8.5e4, 1.5e5),
  ssc_sd = runif(1, 0.5, 1.0),
  fl1_mean = runif(1, 1e5, 1e6),
  fl1_sd = runif(1, 0.5, 1.0),
  fl2_mean = runif(1, 1e5, 1e6),
  fl2_sd = runif(1, 0.5, 1.0),
  fl3_mean = runif(1, 1e5, 1e6),
  fl3_sd = runif(1, 0.5, 1.0)
) {
  # Calculate total events per cell definition
  n_live <- round(n_events * live_pct)
  n_debris <- round(n_events * debris_pct)
  n_dead <- n_events - n_live - n_debris
  n_size <- round(n_live * size_pct)
  n_singlets <- round(n_size * singlet_pct)
  n_doublets <- n_size - n_singlets
  n_fl2_positive <- round(n_singlets * fl2_positive_pct)
  n_fl2_negative <- round(n_singlets - n_fl2_positive)
  n_fl3_positive <- round(n_singlets * fl3_positive_pct)
  n_fl3_negative <- round(n_singlets - n_fl3_positive)

  # Initiate data frame with population identifier: live, debris, or dead
  df <- tibble(
    population = c(
      rep("live", max(0, n_live)),
      rep("debris", max(0, n_debris)),
      rep("dead", max(0, n_dead))
    )
  ) |>

    # Create distribution of FSC-A based on cell population
    mutate(
      `FSC-A` = case_when(
        population == "debris" ~ rnorm(n(), 1e4, 1e4),
        .default = rgamma(n(), shape = 4, rate = 4 / fsc_mean)
      )
    ) |>

    # Create distribution of FL1-A (live/dead dye) based on cell population
    mutate(
      `FL1-A` = case_when(
        population == "live" ~ rlnorm(n(), log(5e4), 0.5),
        population == "debris" ~ rlnorm(n(), log(3.5e4), 0.5),
        population == "dead" ~ rlnorm(n(), log(fl1_mean), fl1_sd),
        .default = NA_real_
      )
    ) |>

    # Create distribution of SSC-A based on cell population
    mutate(
      `SSC-A` = case_when(
        population == "debris" ~ rnorm(n(), 1e4, 1e4),
        .default = rgamma(n(), shape = 4, rate = 4 / ssc_mean) +
          (1e-1 *
            `FSC-A`)
      )
    )
  df
}
