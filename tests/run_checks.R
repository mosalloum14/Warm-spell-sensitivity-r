# From the repository root: Rscript tests/run_checks.R
# Recomputes outputs, then tests event accounting, thresholds, HAC and file hashes.
if (!requireNamespace("sandwich", quietly = TRUE))
  stop("Install dependencies first: source('install_dependencies.R').")
source("run_analysis.R")
same <- function(a, b) stopifnot(isTRUE(all.equal(a, b, tolerance = 1e-10,
                                                check.attributes = FALSE)))
checks <- character()
record <- function(label) { checks <<- c(checks, label); message("PASS: ", label) }

# Qualify a cross-year run before allocating its days to calendar years.
event_fixture <- function() {
  saved_years <- years
  on.exit(assign("years", saved_years, envir = .GlobalEnv))
  assign("years", 2000:2001, envir = .GlobalEnv)
  dates <- seq(as.Date("2000-12-28"), as.Date("2001-01-05"), by = "day")
  ev <- find_events(c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))
  a <- annual_metrics(ev, dates, 4L)
  same(a$n_events, c(1, 0)); same(a$event_days, c(2, 2))
  same(a$mean_duration[1], 4); stopifnot(is.na(a$mean_duration[2]))
  none <- annual_metrics(ev, dates, 5L)
  same(none$n_events, c(0, 0)); same(none$event_days, c(0, 0))
  stopifnot(all(is.na(none$mean_duration)))
  edge <- find_events(c(TRUE, TRUE, rep(FALSE, 5), TRUE, TRUE))
  cut <- annual_metrics(edge, dates, 4L)
  stopifnot(all(is.na(cut$n_events)), all(is.na(cut$event_days)))
}
event_fixture()
record("Cross-year counts/days, zero-event duration and record-edge exclusions")

# A leap-day temperature must not enter the reference quantile sample.
fixture <- data.frame(date = expected_dates,
                      tmax = as.numeric(format(expected_dates, "%j")))
q0 <- daily_thresholds(fixture)
leap <- format(fixture$date, "%m-%d") == "02-29"
fixture$tmax[leap] <- 1e6
same(q0, daily_thresholds(fixture))
if (any(leap)) {
  i <- which(leap)[1]
  same(q0[i, ], (q0[i - 1L, ] + q0[i + 1L, ]) / 2)
}
record("Leap-day exclusion and interpolated threshold")

# Compare complete annual series with an external covariance implementation.
# Gapped years use a separate full 2x2 covariance calculation, with explicit
# calendar-lag pairs rather than compressing gaps or using the scalar formula.
calendar_covariance <- function(x, value, lag) {
  design <- cbind(1, x - mean(x))
  fit <- lm.fit(design, value)
  scores <- design * as.numeric(fit$residuals)
  meat <- crossprod(scores)
  if (lag > 0L) for (j in seq_len(lag)) {
    pairs <- which(outer(x, x, "-") == j, arr.ind = TRUE)
    if (nrow(pairs)) {
      cross <- crossprod(scores[pairs[, 1], , drop = FALSE],
                         scores[pairs[, 2], , drop = FALSE])
      meat <- meat + (1 - j / (lag + 1)) * (cross + t(cross))
    }
  }
  bread <- solve(crossprod(design))
  n <- length(value)
  (bread %*% meat %*% bread) * n / (n - 2)
}
complete_checks <- gap_checks <- 0L
groups <- split(annual_table,
  annual_table[c("lon", "lat", "percentile", "min_days")], drop = TRUE)
for (a in groups) {
  a <- a[match(years, a$year), ]
  for (metric in metrics) {
    v <- a[[metric]]; ok <- is.finite(v)
    result <- fit_trend(v, lag_base)
    if (!is.finite(result$se_per_decade)) next
    x <- years[ok]; value <- v[ok]
    if (all(ok)) {
      model <- lm(value ~ I(x - mean(x)))
      covariance <- sandwich::NeweyWest(model, lag = lag_base,
                                        prewhite = FALSE, adjust = TRUE)
      complete_checks <- complete_checks + 1L
    } else {
      covariance <- calendar_covariance(x, value, lag_base)
      gap_checks <- gap_checks + 1L
    }
    se <- 10 * sqrt(covariance[2, 2])
    same(result$se_per_decade, se)
    margin <- qt((1 + confidence_level) / 2, sum(ok) - 2) * se
    same(c(result$ci_low, result$ci_high),
         result$change_per_decade + c(-margin, margin))
    same(result$p_value,
         2 * pt(-abs(result$change_per_decade / se), sum(ok) - 2))
  }
}
stopifnot(complete_checks > 0L, gap_checks > 0L)
record(sprintf("HAC: %d complete series vs sandwich; %d gapped series vs calendar covariance",
               complete_checks, gap_checks))

# Check the minimum-year guard without changing the analysis settings.
sparse <- rep(NA_real_, length(years)); sparse[seq_len(min_years - 1L)] <-
  seq_len(min_years - 1L)
guard <- fit_trend(sparse, lag_base)
stopifnot(guard$n_years == min_years - 1L, is.na(guard$change_per_decade))
record("Minimum-year eligibility guard")

# A genuine hash failure must stop input use; test with a disposable copy.
hash_fixture <- function() {
  root <- tempfile("snapshot-check-"); dir.create(file.path(root, "data"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE))
  manifest <- read.csv("data/snapshot_manifest.csv", stringsAsFactors = FALSE)[1, ]
  origin <- manifest$path
  manifest$path <- "data/example.json.gz"
  file.copy(origin, file.path(root, manifest$path))
  write.csv(manifest, file.path(root, "data/snapshot_manifest.csv"), row.names = FALSE)
  same(verify_power_snapshot(root = root, quiet = TRUE), 1L)
  con <- file(file.path(root, manifest$path), "ab")
  writeBin(as.raw(0L), con); close(con)
  failed <- tryCatch({ verify_power_snapshot(root = root, quiet = TRUE); FALSE },
                     error = function(e) grepl("SHA-256 mismatch", conditionMessage(e)))
  stopifnot(failed)
}
hash_fixture()
record("Snapshot digest accepts original bytes and rejects modified bytes")
writeLines(c(sprintf("R %s | sandwich %s | digest %s", getRversion(),
                     packageVersion("sandwich"), packageVersion("digest")), checks),
           "results/checks.txt")
message("All reproducible checks passed. Report: results/checks.txt")
