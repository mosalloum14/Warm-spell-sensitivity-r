library(ggplot2)
options(timeout = max(120, getOption("timeout")))

# 1. Settings
area <- c(xmin = 36, ymin = 32.5, xmax = 42, ymax = 37)
years <- 1981:2025
ref_years <- 1986:2005
percentiles <- c(P90 = 0.90, P95 = 0.95)
min_days <- c(2L, 3L, 4L)
half_window <- 2L
confidence_level <- 0.95
alpha <- 0.05
min_years <- 15L
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
lag_base <- min(length(years) - 2L, floor(4 * (length(years) / 100)^(2 / 9)))
lag_checks <- sort(unique(pmin(length(years) - 2L,
                              round(lag_base * c(0.5, 1, 2)))))

# 2. Data
cache_dir <- file.path("power_cache", paste(c("T2M_MAX", unname(area),
  paste0(min(years), "0101"), paste0(max(years), "1231"), "UTC"), collapse = "_"))
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
snapshot_dir <- file.path("data", "power_snapshot", basename(cache_dir))
snapshot_verified <- FALSE
if (dir.exists(snapshot_dir) && isTRUE(getOption("warmspell.verify_snapshot", TRUE))) {
  source("scripts/verify_snapshot.R")
  verify_power_snapshot(snapshot_dir = snapshot_dir)
  snapshot_verified <- TRUE
}
read_year <- function(y) {
  path <- file.path(cache_dir, paste0(y, ".rds"))
  snapshot <- file.path("data", "power_snapshot", basename(cache_dir), paste0(y, ".json.gz"))
  if (file.exists(snapshot)) {
    con <- gzfile(snapshot, "rb")
    on.exit(close(con), add = TRUE)
    obj <- jsonlite::fromJSON(con, simplifyVector = FALSE)
  } else if (file.exists(path)) obj <- readRDS(path) else {
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

# 3. Thresholds
calendar <- format(seq(as.Date("2001-01-01"), by = "day", length.out = 365), "%m-%d")
daily_thresholds <- function(d, reference = ref_years) {
  md <- format(d$date, "%m-%d")
  ref <- as.integer(format(d$date, "%Y")) %in% reference & md != "02-29"
  by_day <- split(d$tmax[ref], factor(md[ref], levels = calendar))
  q <- bind(lapply(seq_along(calendar), function(i) {
    k <- (i - 1L + seq.int(-half_window, half_window)) %% 365L + 1L
    quantile(unlist(by_day[k]), percentiles, type = 8, names = FALSE)
  }))
  dimnames(q) <- list(calendar, names(percentiles))
  q <- rbind(q, "02-29" = colMeans(q[c("02-28", "03-01"), , drop = FALSE]))
  q[md, , drop = FALSE]
}

# 4. Events
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
  out[years %in% yy[ev$first[ev$cut]], c("n_events", "mean_duration")] <- NA
  out$event_days[years %in% yy[indices(ev[ev$cut, ])]] <- NA
  out
}

# 5. Trend function
fit_trend <- function(value, lag, minimum_years = min_years) {
  ok <- is.finite(value)
  x <- years[ok]
  v <- value[ok]
  n <- length(v)
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

# 6. Annual indices
annual_indices <- function(reference = ref_years) {
  annual <- list()
  for (d in cell_series) {
    d <- d[order(d$date), ]
    stopifnot(identical(as.numeric(d$date), as.numeric(expected_dates)))
    q <- daily_thresholds(d, reference)
    for (p in names(percentiles)) {
      ev <- find_events(d$tmax > q[, p])
      for (m in min_days) {
        info <- data.frame(lon = d$lon[1], lat = d$lat[1], percentile = p, min_days = m)
        annual[[length(annual) + 1L]] <- cbind(info, annual_metrics(ev, d$date, m))
      }
    }
  }
  bind(annual)
}
annual_table <- annual_indices()

# 7. Trend estimates
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
  snapshot_SHA256_verified = snapshot_verified,
  input_priority = "bundled snapshot, then local RDS cache, then API",
  planned_tests = nrow(trends), estimated_tests = sum(is.finite(trends$p_value))),
  "results/analysis_settings.rds")
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")

# 8. Maps
prepare_map <- function(d) {
  d$percentile <- factor(d$percentile, levels = names(percentiles))
  d$min_days <- factor(d$min_days, levels = min_days, labels = paste0(min_days, "+ days"))
  d
}
map_coordinates <- function() coord_fixed(
  ratio = 1 / cos(mean(area[c("ymin", "ymax")]) * pi / 180),
  xlim = range(daily$lon) + c(-0.3125, 0.3125),
  ylim = range(daily$lat) + c(-0.25, 0.25), expand = FALSE)
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
      subtitle = sprintf("%d-%d | reference %d-%d | OLS | BY p < %.2f | HAC lag %d",
                         min(years), max(years), min(ref_years), max(ref_years), alpha, lag_base),
      caption = paste(note, sprintf("Grey: unavailable slope; at least %d years required.", min_years),
                      "Full native cells, including their outer half-cell edges."),
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

# 9. Coverage
coverage <- prepare_map(trends[trends$metric == "mean_duration", ])
p_coverage <- ggplot(coverage, aes(lon, lat, fill = n_years)) +
  geom_tile() +
  geom_point(data = coverage[coverage$n_years < min_years, ], shape = 4, size = 1.5) +
  facet_grid(percentile ~ min_days) +
  scale_fill_viridis_c(limits = c(0, length(years)), name = "Years") +
  scale_x_continuous(n.breaks = 4) + map_coordinates() +
  labs(title = "Years available for conditional mean-duration trends",
    subtitle = sprintf("%d-%d | reference %d-%d | Cross: fewer than %d years",
                       min(years), max(years), min(ref_years), max(ref_years), min_years),
    x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), panel.spacing.x = grid::unit(1, "lines"))
ggsave("results/mean_duration_coverage.png", p_coverage,
       width = 11, height = 7, dpi = 300, bg = "white")
if (interactive()) print(p_coverage)
message("Done: ", length(cell_series), " cells; ", nrow(annual_table),
        " annual rows. Tables and maps are in results/.")
