# A limited cross-product diagnostic, using archived ERA5-Land via Open-Meteo.
# Site selection precedes looking at trends: closest POWER center to the rectangle
# midpoint, cosine-scaled distance, with the southern latitude breaking a tie.
# Different cell footprints prevent interpreting agreement as local validation.
compare_era5_land <- function() {
  folder <- file.path("data", "era5_land_snapshot")
  meta <- jsonlite::fromJSON(file.path(folder, "metadata.json"))
  path <- file.path(folder, "center_cell_1981_2025.json")
  stopifnot(identical(digest::digest(file = path, algo = "sha256"), meta$sha256))
  locations <- unique(daily[c("lon", "lat")])
  locations <- locations[order(locations$lat, locations$lon), ]
  center_lon <- mean(area[c("xmin", "xmax")])
  center_lat <- mean(area[c("ymin", "ymax")])
  distance <- ((locations$lon - center_lon) * cos(center_lat * pi / 180))^2 +
              (locations$lat - center_lat)^2
  site <- locations[which.min(distance), ]
  if (abs(site$lon - meta$requested_longitude) > 1e-8 ||
      abs(site$lat - meta$requested_latitude) > 1e-8 ||
      !"P90" %in% names(percentiles) || percentiles["P90"] != 0.90) {
    message("ERA5-Land diagnostic skipped: settings do not match its archived site/P90.")
    return(invisible(FALSE))
  }
  raw <- jsonlite::fromJSON(path)
  stopifnot(raw$utc_offset_seconds == 0,
            identical(raw$daily_units$temperature_2m_max, "\u00b0C"))
  era <- data.frame(date = as.Date(raw$daily$time),
                    tmax = raw$daily$temperature_2m_max)
  if (!all(expected_dates %in% era$date)) {
    message("ERA5-Land diagnostic skipped: archived dates do not cover this run.")
    return(invisible(FALSE))
  }
  era <- era[match(expected_dates, era$date), ]
  power <- daily[daily$lon == site$lon & daily$lat == site$lat, c("date", "tmax")]
  power <- power[match(expected_dates, power$date), ]
  stopifnot(!anyNA(era), !anyNA(power),
            identical(as.numeric(era$date), as.numeric(expected_dates)),
            identical(as.numeric(power$date), as.numeric(expected_dates)))
  annual_comparison <- bind(lapply(c("POWER", "ERA5-Land"), function(product) {
    d <- if (product == "POWER") power else era
    yy <- factor(as.integer(format(d$date, "%Y")), levels = years)
    q <- daily_thresholds(d)
    spells <- annual_metrics(find_events(d$tmax > q[, "P90"]), d$date, 3L)
    data.frame(product = product, year = years,
      mean_tmax = as.vector(tapply(d$tmax, yy, mean)),
      n_events_P90_3days = spells$n_events,
      requested_lon = site$lon, requested_lat = site$lat,
      grid_lon = if (product == "POWER") site$lon else raw$longitude,
      grid_lat = if (product == "POWER") site$lat else raw$latitude)
  }))
  diagnostic_metrics <- c("mean_tmax", "n_events_P90_3days")
  diagnostic_trends <- bind(lapply(split(annual_comparison,
    annual_comparison$product), function(d) bind(lapply(diagnostic_metrics, function(metric) {
      fit <- fit_trend(d[[metric]], lag_base)
      data.frame(product = d$product[1], metric = metric,
        fit[c("n_years", "change_per_decade", "sen_per_decade")])
    }))))
  write.csv(annual_comparison, "results/cross_product_annual.csv", row.names = FALSE, na = "")
  write.csv(diagnostic_trends, "results/cross_product_trends.csv", row.names = FALSE, na = "")
  plot_data <- bind(lapply(diagnostic_metrics, function(metric)
    data.frame(annual_comparison[c("product", "year")], metric = metric,
               value = annual_comparison[[metric]])))
  plot_data$metric <- factor(plot_data$metric, levels = diagnostic_metrics,
    labels = c("Annual mean daily Tmax (C)", "P90, 3+ days: annual events"))
  p <- ggplot(plot_data, aes(year, value, color = product)) + geom_line(alpha = 0.7) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    labs(title = "One-site POWER / ERA5-Land diagnostic",
      subtitle = sprintf("Requested POWER center: %.2f E, %.2f N | %d-%d, UTC",
                         site$lon, site$lat, min(years), max(years)),
      caption = "Nearest cells have different footprints; lines are OLS fits, not accuracy validation.",
      x = "Year", y = NULL, color = NULL) + theme_minimal(base_size = 11)
  ggsave("results/cross_product_comparison.png", p, width = 9, height = 7,
         dpi = 300, bg = "white")
  if (interactive()) print(p)
  message("One-site ERA5-Land diagnostic saved.")
  invisible(diagnostic_trends)
}
compare_era5_land()
