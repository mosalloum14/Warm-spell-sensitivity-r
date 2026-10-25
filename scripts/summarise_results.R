# Small descriptive tables used in the README; no new statistical tests.
definition_summary <- bind(lapply(split(trends,
  trends[c("percentile", "min_days", "metric")], drop = TRUE), function(d) {
  eligible <- is.finite(d$change_per_decade)
  data.frame(d[1, c("percentile", "min_days", "metric")],
    planned_cells = nrow(d), eligible_cells = sum(eligible),
    median_OLS = median(d$change_per_decade[eligible]),
    median_Sen = median(d$sen_per_decade[is.finite(d$sen_per_decade)]),
    significant_BY = sum(d$p_adjusted_BY < alpha, na.rm = TRUE))
}))
write.csv(definition_summary, "results/definition_summary.csv", row.names = FALSE, na = "")

# Recompute BY independently for each bandwidth, retaining the planned family.
# Restrict to available-year fits so the sampling design is held fixed.
fit_keys <- c("lon", "lat", "percentile", "min_days", "metric")
key <- function(d) do.call(paste, c(d[fit_keys], sep = ":"))
primary_key <- key(trends)
hac_lag_summary <- bind(lapply(lag_checks, function(bandwidth) {
  d <- all_trends[all_trends$sample == "available" &
                  all_trends$lag_years == bandwidth, ]
  d$p_adjusted_BY <- p.adjust(d$p_value, method = "BY", n = nrow(d))
  base <- trends[match(key(d), primary_key), ]
  stopifnot(!anyNA(base$metric))
  significant <- is.finite(d$p_adjusted_BY) & d$p_adjusted_BY < alpha
  base_significant <- is.finite(base$p_adjusted_BY) & base$p_adjusted_BY < alpha
  bind(lapply(metrics, function(metric) {
    take <- d$metric == metric
    data.frame(lag_years = bandwidth, metric = metric,
      planned_cases = sum(take), estimated_tests = sum(is.finite(d$p_value[take])),
      significant_BY = sum(significant[take]),
      decisions_changed_vs_primary = sum(significant[take] != base_significant[take]))
  }))
}))
write.csv(hac_lag_summary, "results/hac_lag_sensitivity.csv", row.names = FALSE, na = "")
message("Definition and HAC-bandwidth summaries saved.")
