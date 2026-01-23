plot_experiment <- function(experiment) {
  experiment |>
    bind_rows(.id = "sample") |>
    ggplot(aes(`FSC-A`, `FL1-A`)) +
    geom_hex(bins = 250) +
    facet_grid(Live ~ sample) +
    scale_y_continuous(transform = "asinh")

  experiment |>
    bind_rows(.id = "sample") |>
    filter(Live) |>
    ggplot(aes(`FSC-A`, `SSC-A`)) +
    geom_hex(bins = 250) +
    facet_grid(Size ~ sample)

  experiment |>
    bind_rows(.id = "sample") |>
    filter(Size) |>
    ggplot(aes(`FSC-A`, `FSC-H`)) +
    geom_hex(bins = 250) +
    facet_grid(Singlet ~ sample)

  experiment |>
    bind_rows(.id = "sample") |>
    filter(Singlet) |>
    pivot_longer(contains("positive")) |>
    ggplot(aes(`FL2-A`, `FL3-A`)) +
    geom_hex(bins = 250) +
    facet_wrap(~sample) +
    scale_y_continuous(transform = "asinh") +
    scale_x_continuous(transform = "asinh")
}
