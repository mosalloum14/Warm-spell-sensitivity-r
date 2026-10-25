# Fixed-reference warm spells: an extension of the MSc methodology.
# Full-year local anomalies, not standard WSDI or a health-risk threshold.
# OLS/HAC inference is approximate, conditional on this product and fixed thresholds.
# Internal robustness checks are not validation against independent observations.
# First run only: install.packages(c("jsonlite", "ggplot2"))
library(ggplot2)
options(timeout = max(120, getOption("timeout")))

# 1. Settings: study rectangle and six definitions retained from the MSc.
area <- c(xmin = 36, ymin = 32.5, xmax = 42, ymax = 37)
years <- 1981:2025
ref_years <- 1991:2020
percentiles <- c(P90 = 0.90, P95 = 0.95)
min_days <- c(2L, 3L, 4L)
half_window <- 2L
confidence_level <- 0.95
alpha <- 0.05                          # FDR decision threshold
min_years <- 15L                       # Pragmatic guard against sparse annual samples;
                                      # not a universal guarantee of HAC adequacy.
metrics <- c("n_events", "event_days", "mean_duration")
bind <- function(x) do.call(rbind, x)
expected_dates <- seq(as.Date(paste0(min(years), "-01-01")),
                      as.Date(paste0(max(years), "-12-31")), by = "day")
stopifnot(length(years) >= 3, all(diff(years) == 1),
          length(ref_years) > 0, all(ref_years %in% years),
          all(percentiles > 0 & percentiles < 1),
          !is.null(names(percentiles)), !anyDuplicated(names(percentiles)),
          all(min_days >= 1 & min_days == as.integer(min_days)),
          !anyDuplicated(min_days), half_window >= 0,
          half_window == as.integer(half_window), half_window < 183,
          min_years >= 3, min_years == as.integer(min_years),
          confidence_level > 0, confidence_level < 1, alpha > 0, alpha < 1)
# Prespecified annual-lag rule; half/double bandwidths are sensitivity checks.
lag_base <- min(length(years) - 2L, floor(4 * (length(years) / 100)^(2 / 9)))
lag_checks <- sort(unique(pmin(length(years) - 2L,
                              round(lag_base * c(0.5, 1, 2)))))

# 2. Read cached POWER data, the bundled snapshot, or download from the API.
cache_dir <- file.path("power_cache", paste(c("T2M_MAX", unname(area),
  paste0(min(years), "0101"), paste0(max(years), "1231"), "UTC"), collapse = "_"))
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
read_year <- function(y) {
  path <- file.path(cache_dir, paste0(y, ".rds"))
  snapshot <- file.path("data", "power_snapshot", basename(cache_dir), paste0(y, ".json.gz"))
  if (file.exists(path)) obj <- readRDS(path) else if (file.exists(snapshot)) {
    con <- gzfile(snapshot, "rb")
    on.exit(close(con), add = TRUE)
    obj <- jsonlite::fromJSON(con, simplifyVector = FALSE)
  } else {
    url <- sprintf(paste0(
      "https://power.larc.nasa.gov/api/temporal/daily/regional?",
      "parameters=T2M_MAX&community=AG&format=JSON&time-standard=UTC",
      "&longitude-min=%s&longitude-max=%s&latitude-min=%s&latitude-max=%s",
      "&start=%d0101&end=%d1231"),
      area["xmin"], area["xmax"], area["ymin"], area["ymax"], y, y)
    obj <- jsonlite::fromJSON(url, simplifyVector = FALSE)
  }
  stopifnot(identical(obj$header$time_standard, "UTC"),
            identical(obj$parameters$T2M_MAX$units, "C"), length(obj$features) > 0)
  if (!file.exists(path)) saveRDS(obj, path)
  bind(lapply(obj$features, function(f) {
    v <- unlist(f$properties$parameter$T2M_MAX)
    v[v == obj$header$fill_value] <- NA_real_
    data.frame(lon = f$geometry$coordinates[[1]], lat = f$geometry$coordinates[[2]],
               date = as.Date(names(v), "%Y%m%d"), tmax = as.numeric(v))
  }))
}
daily <- bind(lapply(years, function(y) { message("Reading ", y); read_year(y) }))
stopifnot(!anyNA(daily), all(is.finite(daily$tmax)))
cell_series <- split(daily, daily[c("lon", "lat")], drop = TRUE)

# 3. Calendar-day quantiles: fixed reference, circular five-day window, type 8.
calendar <- format(seq(as.Date("2001-01-01"), by = "day", length.out = 365), "%m-%d")
daily_thresholds <- function(d) {
  md <- format(d$date, "%m-%d")
  ref <- as.integer(format(d$date, "%Y")) %in% ref_years & md != "02-29"
  by_day <- split(d$tmax[ref], factor(md[ref], levels = calendar))
  q <- bind(lapply(seq_along(calendar), function(i) {
    k <- (i - 1L + seq.int(-half_window, half_window)) %% 365L + 1L
    quantile(unlist(by_day[k]), percentiles, type = 8, names = FALSE)
  }))
  dimnames(q) <- list(calendar, names(percentiles))
  q <- rbind(q, "02-29" = colMeans(q[c("02-28", "03-01"), , drop = FALSE]))
  q[md, , drop = FALSE]
}

# 4. Detect complete runs before grouping years; do not split at New Year.
find_events <- function(hot) {
  r <- rle(hot)
  last <- cumsum(r$lengths)
  first <- last - r$lengths + 1L
  data.frame(first = first[r$values], last = last[r$values],
             days = r$lengths[r$values],
             cut = (first == 1L | last == length(hot))[r$values])
}
annual_metrics <- function(ev, date, min_d) {
  yy <- as.integer(format(date, "%Y"))
  e <- ev[ev$days >= min_d & !ev$cut, ]
  onset <- factor(yy[e$first], levels = years)
  indices <- function(e) unlist(Map(seq.int, e$first, e$last), use.names = FALSE)
  out <- data.frame(year = years,
    n_events = as.vector(table(onset)),
    event_days = tabulate(match(yy[indices(e)], years), nbins = length(years)),
    mean_duration = as.vector(tapply(e$days, onset, mean)))
  # Count/full duration belong to onset year; event_days belong to actual day-year.
  # Conservatively exclude years affected by any record-edge run, even a short one.
  out[years %in% yy[ev$first[ev$cut]], c("n_events", "mean_duration")] <- NA
  out$event_days[years %in% yy[indices(ev[ev$cut, ])]] <- NA
  out
}

# 5. OLS change/decade, Sen robustness estimate, approximate HAC interval.
# Bartlett weights use actual year distances: missing event-years are not compressed.
fit_trend <- function(value, lag, minimum_years = min_years) {
  ok <- is.finite(value)
  x <- years[ok]; v <- value[ok]; n <- length(v)
  out <- data.frame(n_years = n,
    first_year = if (n) min(x) else NA_integer_,
    last_year = if (n) max(x) else NA_integer_,
    change_per_decade = NA_real_, sen_per_decade = NA_real_,
    se_per_decade = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_)
  if (n < minimum_years) return(out)
  pair_slopes <- outer(v, v, "-") / outer(x, x, "-")
  out$sen_per_decade <- 10 * median(pair_slopes[upper.tri(pair_slopes)])
  xc <- x - mean(x)
  fit <- lm(v ~ xc)
  out$change_per_decade <- if (length(unique(v)) == 1L) 0 else 10 * unname(coef(fit)[2])
  if (max(abs(residuals(fit))) <=
      sqrt(.Machine$double.eps) * max(1, max(abs(v)))) return(out)
  score <- xc * residuals(fit) / sum(xc^2)
  kernel <- pmax(1 - abs(outer(x, x, "-")) / (lag + 1), 0)
  se <- 10 * sqrt(max(0, sum(score * (kernel %*% score)) * n / (n - 2)))
  if (se <= 0) return(out)
  out$se_per_decade <- se
  margin <- qt((1 + confidence_level) / 2, df = n - 2) * se
  out$ci_low <- out$change_per_decade - margin
  out$ci_high <- out$change_per_decade + margin
  out$p_value <- 2 * pt(-abs(out$change_per_decade / se), df = n - 2)
  out
}

# 6. Annual indices for every cell and definition; exact dates are mandatory.
annual <- list()
for (d in cell_series) {
  d <- d[order(d$date), ]
  stopifnot(identical(as.numeric(d$date), as.numeric(expected_dates)))
  q <- daily_thresholds(d)
  for (p in names(percentiles)) {
    ev <- find_events(d$tmax > q[, p])
    for (m in min_days) {
      info <- data.frame(lon = d$lon[1], lat = d$lat[1], percentile = p, min_days = m)
      annual[[length(annual) + 1L]] <- cbind(info, annual_metrics(ev, d$date, m))
    }
  }
}
annual_table <- bind(annual)

# 7. Available years are primary. Common event-years are a secondary duration check.
# Mean duration is conditional on events; a no-event year has NA duration, not zero.
fits <- list()
for (cell in split(annual_table, annual_table[c("lon", "lat")], drop = TRUE)) {
  definitions <- split(cell, cell[c("percentile", "min_days")], drop = TRUE)
  common <- Reduce(intersect, lapply(definitions, function(a)
    a$year[is.finite(a$mean_duration)]))
  for (a in definitions) {
    a <- a[match(years, a$year), ]
    info <- a[1, c("lon", "lat", "percentile", "min_days")]
    for (metric in metrics) for (lag in lag_checks) {
      fits[[length(fits) + 1L]] <- data.frame(info, metric = metric,
        sample = "available", lag_years = lag, fit_trend(a[[metric]], lag))
    }
    v <- a$mean_duration
    v[!years %in% common] <- NA
    fits[[length(fits) + 1L]] <- data.frame(info, metric = "mean_duration",
      sample = "common", lag_years = lag_base, fit_trend(v, lag_base))
  }
}
all_trends <- bind(fits)
trends <- subset(all_trends, sample == "available" & lag_years == lag_base)
# One prespecified family: all cells x six definitions x three primary metrics.
# Unestimable cases remain in the planned family; report valid p-values separately.
# BY accommodates dependence; inference remains approximate.
trends$p_adjusted_BY <- p.adjust(trends$p_value, method = "BY", n = nrow(trends))
trend_checks <- subset(all_trends, sample == "common" | lag_years != lag_base)
write.csv(annual_table, "results/heatwave_annual.csv", row.names = FALSE, na = "")
write.csv(trends, "results/heatwave_trends.csv", row.names = FALSE, na = "")
write.csv(trend_checks, "results/heatwave_trend_checks.csv", row.names = FALSE, na = "")
saveRDS(list(area = area, years = years, ref_years = ref_years,
  percentiles = percentiles, min_days = min_days, half_window = half_window,
  min_years = min_years, confidence_level = confidence_level, FDR_level = alpha,
  minimum_years_rule = "15-year default is a pragmatic sparse-sample guard, not a universal standard",
  time_standard = "UTC", variable = "T2M_MAX", units = "C",
  primary_HAC_lag = lag_base, HAC_lag_checks = lag_checks,
  inference = "Approximate OLS/HAC; conditional on fixed thresholds and POWER data",
  multiple_testing = "BY across all planned primary cases, including unestimable cases",
  planned_tests = nrow(trends), estimated_tests = sum(is.finite(trends$p_value))),
  "results/analysis_settings.rds")
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")

# 8. Identical scales across definitions within each metric; dots use primary BY p.
prepare_map <- function(d) {
  d$percentile <- factor(d$percentile, levels = names(percentiles))
  d$min_days <- factor(d$min_days, levels = min_days, labels = paste0(min_days, "+ days"))
  d
}
map_coordinates <- function() coord_fixed(
  ratio = 1 / cos(mean(area[c("ymin", "ymax")]) * pi / 180),
  xlim = unname(area[c("xmin", "xmax")]),
  ylim = unname(area[c("ymin", "ymax")]), expand = FALSE)
plot_trend <- function(metric, title, unit, note) {
  d <- prepare_map(trends[trends$metric == metric, ])
  z <- abs(d$change_per_decade[is.finite(d$change_per_decade)])
  lim <- if (length(z) && max(z) > 0) max(z) else 1
  p <- ggplot(d, aes(lon, lat, fill = change_per_decade)) +
    geom_tile() +
    geom_point(data = d[which(d$p_adjusted_BY < alpha), ], size = 0.55) +
    facet_grid(percentile ~ min_days) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      limits = c(-lim, lim), na.value = "grey80", name = unit) +
    scale_x_continuous(n.breaks = 4) +
    map_coordinates() +
    labs(title = title,
      subtitle = sprintf("%d-%d | OLS | dots: global BY-adjusted p < %.2f (HAC lag %d)",
                         min(years), max(years), alpha, lag_base),
      caption = paste(note, sprintf("Grey: unavailable slope; at least %d years required.", min_years),
                      "Cells clipped to the study rectangle."),
      x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), panel.spacing.x = grid::unit(1, "lines"))
  ggsave(file.path("results", paste0(metric, "_trends.png")), p,
         width = 11, height = 7, dpi = 300, bg = "white")
  if (interactive()) print(p)
}
plot_trend("n_events", "Change in warm-spell frequency",
           "Events/year\nchange per decade", "Events are counted by onset year.")
plot_trend("event_days", "Change in annual days in warm spells",
           "Days/year\nchange per decade", "Each qualifying day belongs to its calendar year.")
plot_trend("mean_duration", "Change in mean warm-spell duration",
           "Days/event\nchange per decade", "Conditional on events; event-bearing years vary by definition.")

# 9. Available event-bearing years, including those excluded by min_years.
coverage <- prepare_map(trends[trends$metric == "mean_duration", ])
p_coverage <- ggplot(coverage, aes(lon, lat, fill = n_years)) +
  geom_tile() +
  geom_point(data = coverage[coverage$n_years < min_years, ], shape = 4, size = 1.5) +
  facet_grid(percentile ~ min_days) +
  scale_fill_viridis_c(limits = c(0, length(years)), name = "Years") +
  scale_x_continuous(n.breaks = 4) + map_coordinates() +
  labs(title = "Years available for conditional mean-duration trends",
    subtitle = sprintf("%d-%d | Cross: fewer than %d years; trend omitted",
                       min(years), max(years), min_years),
    x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), panel.spacing.x = grid::unit(1, "lines"))
ggsave("results/mean_duration_coverage.png", p_coverage,
       width = 11, height = 7, dpi = 300, bg = "white")
if (interactive()) print(p_coverage)
message("Done: ", length(cell_series), " cells; ", nrow(annual_table),
        " annual rows. Tables and maps are in results/.")

# Method references:
# https://www.climdex.org/learn/indices/ (standard indices differ from our spell lengths)
# https://search.r-project.org/CRAN/refmans/sandwich/html/NeweyWest.html
# https://stat.ethz.ch/R-manual/R-patched/library/stats/html/p.adjust.html
