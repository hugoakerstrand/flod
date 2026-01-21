#' @export
synthetic_sample <- function(
  n_events = 20000,
  debris_pct = runif(1, 0.01, 0.15),
  live_pct = runif(1, 0.75, 0.99),
  singlet_pct = runif(1, 0.85, 0.99),
  fl2_positive_pct = runif(1, 0.05, 0.95),
  fl3_positive_pct = runif(1, 0.05, 0.95),
  debris_sd = runif(1, 1e4, 1e5),
  fsc_mean = runif(1, 1e5, 2e5),
  ssc_mean = runif(1, 8.5e4, 1.5e5),
  fl1_mean = runif(1, 1e5, 1e6),
  fl1_sd = runif(1, 0.5, 2.0),
  fl2_mean = runif(1, 0.5e5, 1e6),
  fl2_sd = runif(1, 0.1, 2.0),
  fl3_mean = runif(1, 0.5e5, 1e6),
  fl3_sd = runif(1, 0.1, 2.0)
) {
  # Calculate total events per cell definition
  n_live <- round(n_events * live_pct)
  n_debris <- round(n_events * debris_pct)
  n_dead <- n_events - n_live - n_debris

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
        population == "debris" ~ rnorm(n(), 1e4, 2e4),
        .default = rgamma(n(), shape = 4, rate = 4 / fsc_mean)
      )
    ) |>

    # Create distribution of FL1-A (live/dead dye) based on cell population
    mutate(
      `FL1-A` = case_when(
        population == "live" ~ rlnorm(n(), log(5), 0.5),
        population == "debris" ~ rlnorm(n(), log(3), 0.5),
        population == "dead" ~ rlnorm(n(), log(fl1_mean), fl1_sd),
        .default = NA_real_
      )
    ) |>

    # Create distribution of SSC-A based on cell population
    mutate(
      `SSC-A` = case_when(
        population == "debris" ~ rnorm(n(), 1e4, debris_sd),
        .default = rgamma(n(), shape = 4, rate = 4 / ssc_mean) +
          (1e-1 *
            `FSC-A`)
      )
    ) |>

    # Create distribution of FSC-H based on cell population
    mutate(
      # Create singlet indicator
      singlet = runif(n()) < singlet_pct,

      # Generate FSC-H based on singlet status
      `FSC-H` = case_when(
        singlet ~ `FSC-A` *
          rnorm(n(), mean = 1, sd = 0.15 / pmax(`FSC-A` / fsc_mean, 0.1)),
        .default = `FSC-A` *
          rnorm(n(), mean = 0.3 + 0.2 * pmax(`FSC-A` / fsc_mean, 0), sd = 0.1),
      )
    ) |>

    # Create distribution of marker positive and negative given their respective pct and mean
    mutate(
      `FL2_positive` = runif(n()) < fl2_positive_pct,
      `FL3_positive` = runif(n()) < fl3_positive_pct,

      # Generate FL2-A based on FL2_positive status
      `FL2-A` = if_else(
        FL2_positive,
        rlnorm(n(), meanlog = log(fl2_mean), sdlog = fl2_sd),
        rlnorm(n(), meanlog = log(5), sdlog = fl2_sd)
      ),

      # Generate FL3-A based on FL3_positive status
      `FL3-A` = if_else(
        FL3_positive,
        rlnorm(n(), meanlog = log(fl3_mean), sdlog = fl3_sd),
        rlnorm(n(), meanlog = log(5), sdlog = fl3_sd)
      )
    )

  df
}
