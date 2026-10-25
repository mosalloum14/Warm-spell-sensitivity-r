# Period and reference sensitivity
check_period_reference <- function() {
  references <- list(ref_years, 1991:2020)
  label <- function(x) paste(min(x), max(x), sep = "-")
  names(references) <- vapply(references, label, character(1))
  references <- references[!duplicated(names(references))]
  windows <- list(years[years <= 2017L], years)
  if (length(windows[[1]]) < 3L || identical(windows[[1]], years) ||
      !all(unlist(references) %in% years)) {
    message("Period/reference comparison skipped: required years are unavailable.")
    return(invisible(NULL))
  }
  names(windows) <- vapply(windows, label, character(1))
  columns <- c("lon", "lat", "percentile", "min_days")
  estimates <- list()
  for (reference in names(references)) {
    a <- if (identical(references[[reference]], ref_years)) annual_table else
      annual_indices(references[[reference]])
    groups <- split(a, a[columns], drop = TRUE)
    for (window in names(windows)) for (d in groups) {
      d <- d[match(years, d$year), ]
      for (metric in metrics) {
        value <- d[[metric]]
        value[!years %in% windows[[window]]] <- NA_real_
        fit <- fit_trend(value, lag_base)
        estimates[[length(estimates) + 1L]] <- data.frame(
          reference = reference, trend_window = window, d[1, columns],
          metric = metric, fit[c("n_years", "first_year", "last_year",
                                "change_per_decade", "sen_per_decade")])
      }
    }
  }
  estimates <- bind(estimates)
  mid <- function(x) if (any(is.finite(x))) median(x[is.finite(x)]) else NA_real_
  summary <- bind(lapply(split(estimates,
    estimates[c("percentile", "min_days", "metric")], drop = TRUE), function(d) {
    cell <- paste(d$lon, d$lat, sep = ":")
    common <- names(which(tapply(is.finite(d$change_per_decade), cell, all)))
    d$common_cell <- cell %in% common
    bind(lapply(split(d, d[c("reference", "trend_window")], drop = TRUE), function(g) {
      data.frame(g[1, c("reference", "trend_window", "percentile", "min_days", "metric")],
        eligible_cells = sum(is.finite(g$change_per_decade)),
        common_cells = sum(g$common_cell), median_OLS = mid(g$change_per_decade),
        median_Sen = mid(g$sen_per_decade),
        common_median_OLS = mid(g$change_per_decade[g$common_cell]),
        common_median_Sen = mid(g$sen_per_decade[g$common_cell]))
    }))
  }))
  write.csv(estimates, "results/period_reference_trends.csv", row.names = FALSE, na = "")
  write.csv(summary, "results/period_reference_summary.csv", row.names = FALSE, na = "")
  d <- summary
  d$reference <- factor(d$reference, levels = names(references))
  d$trend_window <- factor(d$trend_window, levels = names(windows))
  d$metric <- factor(d$metric, levels = metrics,
    labels = c("Frequency (events/year)", "Event days (days/year)", "Duration (days/event)"))
  p <- ggplot(d, aes(min_days, common_median_OLS, color = reference,
                     linetype = trend_window,
                     group = interaction(reference, trend_window))) +
    geom_line(na.rm = TRUE) + geom_point(na.rm = TRUE) +
    facet_grid(metric ~ percentile, scales = "free_y") +
    scale_x_continuous(breaks = min_days) +
    labs(title = "Trend-window and reference-period sensitivity",
      subtitle = "Same eligible cells across all cases within each definition and metric",
      caption = "Descriptive OLS slopes; sample coverage and Sen slopes are in the CSV tables.",
      x = "Minimum event length (days)", y = "Median change per decade",
      color = "Reference", linetype = "Trend window") +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  ggsave("results/period_reference_comparison.png", p, width = 10, height = 9,
         dpi = 300, bg = "white")
  if (interactive()) print(p)
  message("Period/reference comparison saved.")
  invisible(list(trends = estimates, summary = summary))
}
period_reference_results <- check_period_reference()
