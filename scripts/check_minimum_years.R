# Reuse the main analysis functions and annual data; no extra downloads.
# Recompute BY for each rule while retaining the same planned testing family.
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

# Display the strictest definition, where sparse event-years are most likely.
strictest <- names(percentiles)[which.max(percentiles)]
d <- minimum_years_trends[minimum_years_trends$metric == "mean_duration" &
       minimum_years_trends$percentile == strictest &
       minimum_years_trends$min_days == max(min_days), ]
d$minimum_years <- factor(d$minimum_years, levels = minimum_options,
                          labels = paste0("Minimum ", minimum_options, " years"))
z <- abs(d$change_per_decade[is.finite(d$change_per_decade)])
lim <- if (length(z) && max(z) > 0) max(z) else 1
p <- ggplot(d, aes(lon, lat, fill = change_per_decade)) +
  geom_tile() + geom_point(data = d[which(d$p_adjusted_BY < alpha), ], size = 0.6) +
  facet_wrap(~ minimum_years, nrow = 1) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
    limits = c(-lim, lim), na.value = "grey80", name = "Days/event\nchange per decade") +
  scale_x_continuous(n.breaks = 4) + map_coordinates() +
  labs(title = sprintf("Minimum-year sensitivity: %s, %d+ days", strictest, max(min_days)),
    subtitle = sprintf("Conditional mean duration | grey: excluded | dots: BY-adjusted p < %.2f", alpha),
    caption = sprintf("The planned testing family remains %d for each rule.",
                       nrow(trends)), x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), panel.spacing.x = grid::unit(1, "lines"))
ggsave("results/minimum_years_comparison.png", p, width = 12, height = 4.8,
       dpi = 300, bg = "white")
if (interactive()) print(p)
message("Minimum-year checks complete: 10, 15, 20 years.")
