# Minimum-year sensitivity
minimum_options <- c(10L, 15L, 20L)
group_columns <- c("lon", "lat", "percentile", "min_days")
groups <- split(annual_table, annual_table[group_columns], drop = TRUE)
minimum_years_trends <- bind(lapply(minimum_options, function(cutoff) {
  d <- bind(lapply(groups, function(a) {
    a <- a[match(years, a$year), ]
    info <- a[1, group_columns]
    bind(lapply(metrics, function(metric) data.frame(info, metric = metric,
      minimum_years = cutoff, fit_trend(a[[metric]], lag_base, cutoff))))
  }))
  d$p_adjusted_BY <- p.adjust(d$p_value, method = "BY", n = nrow(d))
  d
}))
minimum_years_summary <- bind(lapply(split(minimum_years_trends,
  minimum_years_trends[c("minimum_years", "metric")], drop = TRUE), function(d) {
  eligible <- is.finite(d$change_per_decade)
  data.frame(d[1, c("minimum_years", "metric")], planned = nrow(d),
    eligible = sum(eligible), estimated_tests = sum(is.finite(d$p_value)),
    positive_OLS = sum(d$change_per_decade > 0, na.rm = TRUE),
    negative_OLS = sum(d$change_per_decade < 0, na.rm = TRUE),
    median_OLS = if (any(eligible)) median(d$change_per_decade[eligible]) else NA_real_,
    significant_BY = sum(d$p_adjusted_BY < alpha, na.rm = TRUE))
}))
write.csv(minimum_years_trends, "results/minimum_years_trends.csv", row.names = FALSE, na = "")
write.csv(minimum_years_summary, "results/minimum_years_sensitivity.csv", row.names = FALSE, na = "")

message("Minimum-year checks complete: 10, 15, 20 years.")
