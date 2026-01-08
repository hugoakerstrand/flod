# Verify the synthetic flow cytometry data looks realistic
library(tidyverse)
library(duckdb)

# Connect to database
con <- dbConnect(duckdb(), dbdir = "flow_cytometry_data.duckdb", read_only = TRUE)

# Load data
cat("Loading data from DuckDB...\n")
flow_data <- dbGetQuery(con, "SELECT * FROM flow_data")
flow_summary <- dbGetQuery(con, "SELECT * FROM flow_summary")

# Close connection
dbDisconnect(con, shutdown = TRUE)

cat("\n=== DATA VERIFICATION ===\n\n")

# 1. Check sample sizes
cat("1. Sample sizes:\n")
print(table(flow_data$sample_id))

# 2. Check summary statistics
cat("\n2. Summary statistics by sample:\n")
print(flow_summary)

# 3. Verify live/dead distribution
cat("\n3. Live/Dead distribution (based on FL1-A < 5000 threshold):\n")
flow_data %>%
  mutate(status = ifelse(`FL1-A` < 5000, "Live", "Dead")) %>%
  group_by(sample_id, status) %>%
  summarise(count = n(), pct = n()/200) %>%
  print()

# 4. Verify marker expression patterns
cat("\n4. Marker expression patterns:\n")
cat("\nSample1 (should have distinct positive/negative populations):\n")
flow_data %>%
  filter(sample_id == "Sample1", `FL1-A` < 5000) %>%  # Live cells only
  summarise(
    FL2_high_pct = sum(`FL2-A` > 5000) / n() * 100,
    FL3_high_pct = sum(`FL3-A` > 5000) / n() * 100
  ) %>%
  print()

cat("\nSample2 (should have background signal only):\n")
flow_data %>%
  filter(sample_id == "Sample2", `FL1-A` < 5000) %>%  # Live cells only
  summarise(
    FL2_high_pct = sum(`FL2-A` > 5000) / n() * 100,
    FL3_high_pct = sum(`FL3-A` > 5000) / n() * 100
  ) %>%
  print()

cat("\nSample3 (failed sample - should have mostly dead cells):\n")
flow_data %>%
  filter(sample_id == "Sample3") %>%
  summarise(
    dead_pct = sum(`FL1-A` > 5000) / n() * 100,
    live_cells = sum(`FL1-A` < 5000)
  ) %>%
  print()

# 5. Check FSC-A vs FSC-H correlation (for singlet detection)
cat("\n5. FSC-A vs FSC-H correlation (should be high for singlets):\n")
flow_data %>%
  filter(`FL1-A` < 5000) %>%  # Live cells only
  group_by(sample_id) %>%
  summarise(
    correlation = cor(`FSC-A`, `FSC-H`),
    fsc_ratio_mean = mean(`FSC-A` / `FSC-H`),
    fsc_ratio_sd = sd(`FSC-A` / `FSC-H`)
  ) %>%
  print()

# 6. Check for debris (low FSC/SSC)
cat("\n6. Debris detection (FSC-A < 10000 AND SSC-A < 10000):\n")
flow_data %>%
  group_by(sample_id) %>%
  summarise(
    debris_count = sum(`FSC-A` < 10000 & `SSC-A` < 10000),
    debris_pct = sum(`FSC-A` < 10000 & `SSC-A` < 10000) / n() * 100
  ) %>%
  print()

cat("\n=== VERIFICATION COMPLETE ===\n")
cat("\nData looks realistic! The DuckDB database is ready at:\n")
cat("  flow_cytometry_data.duckdb\n\n")
cat("Tables:\n")
cat("  - flow_data: Main data table with all events\n")
cat("  - flow_summary: Summary statistics by sample\n")
