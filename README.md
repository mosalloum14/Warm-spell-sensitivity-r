# Warm-spell definition sensitivity in R

An MSc methodology extension by **Mohammad Salloum**: how did the frequency and duration of relative warm spells change during 1981–2025, and how sensitive are the results to the percentile and minimum event length?

The workflow combines open climate data, calendar checks, consecutive-event detection, annual indices, trend estimation and explicit sensitivity analysis. It uses base R, **jsonlite** and **ggplot2**.

## Connection to the MSc thesis

The project extends the percentile-and-persistence framework in *Climate Change Consequences and Scenarios: Syria as Case Study* (University of Tübingen, 2021, section 3.3.5, p. 45). It retains the **P90/P95 thresholds and 2/3/4-day event lengths**, together with the original extraction rectangle.

The data and estimand have changed: the thesis describes daily mean temperature series and a 1986–2005 baseline; this project studies historical **daily maximum temperature**, with a fixed **1991–2020** reference. It is a methodological extension, rather than a numerical reproduction or independent validation of the thesis results.

## Results

The checked run contains **1,643,600 daily records**, **27,000 annual index rows**, and **1,800 planned primary trend cases**. With the 15-year rule, 1,798 primary p-values are estimable.

<!-- RESULTS_TABLE_START -->
| Definition | Frequency OLS | Frequency Sen | Event-day totals OLS | Mean-duration OLS | Duration fits |
| --- | ---: | ---: | ---: | ---: | ---: |
| P90, 2+ days | 2.083 | 2.026 | 8.324 | 0.179 | 100/100 |
| P90, 3+ days | 1.177 | 1.111 | 6.462 | 0.255 | 100/100 |
| P90, 4+ days | 0.779 | 0.702 | 5.254 | 0.294 | 100/100 |
| P95, 2+ days | 1.227 | 1.111 | 4.233 | 0.179 | 100/100 |
| P95, 3+ days | 0.636 | 0.500 | 3.105 | 0.261 | 100/100 |
| P95, 4+ days | 0.420 | 0.000 | 2.536 | 0.323 | 98/100 |

Values are spatial medians of cell-specific changes per decade, not an area-average fitted trend. Metric units are defined below.
<!-- RESULTS_TABLE_END -->

The strongest common finding is an **OLS-based increase in frequency and event-day totals** across all six definitions. Effect size and statistical evidence depend on the definition. In the strictest frequency definition (P95, 4+ days), the spatial median Sen slope is zero despite a positive median OLS slope; agreement between estimators is therefore limited for that sparse count series.

The minimum-year comparison retains 600, 598 and 564 conditional-duration fits under the 10-, 15- and 20-year rules. The coverage loss occurs in P95/4+ days, where 100, 98 and 64 cells remain. Raising the rule from 15 to 20 changes one BY significance decision among shared fits. This check describes sensitivity to sample eligibility, not an independent validation of temperature accuracy.

![Warm-spell frequency trends](results/n_events_trends.png)

![Minimum-year sensitivity](results/minimum_years_comparison.png)

## Conclusion and comparison with the MSc thesis

The results are **qualitatively consistent with the thesis' historical increase in hot days** (1980–2017; pp. 49–51). For 1981–2025, OLS warm-spell frequency slopes are positive in all 600 cell–definition cases; 400 pass the global BY correction. Spatial median changes range from **+0.42 to +2.08 annual events per decade** across definitions. Conditional mean-duration trends are mostly positive, but statistical evidence is less widespread: 60 of 598 eligible cases pass BY.

The thesis describes **daily mean temperature** from [GHCNd station records](https://www.ncei.noaa.gov/products/land-based-station/global-historical-climatology-network-daily) and [Berkeley Earth's station-based gridded reconstruction](https://berkeleyearth.org/methodology/) (pp. 32–33). This project uses **NASA POWER/MERRA-2 daily maximum temperature**, derived from hourly estimates on a 0.5° latitude × 0.625° longitude grid. [Reanalysis combines observations with an atmospheric model](https://power.larc.nasa.gov/docs/methodology/meteorology/); it is not a collection of local station measurements. Differences in variable, reference period (1986–2005 versus 1991–2020), analysis period and spatial treatment prevent direct numerical comparison.

The contribution is a **reproducible, spatially explicit historical extension through 2025**, with calendar-day thresholds and explicit definition and sample-coverage checks. This supports the thesis' general direction of change without verifying its original event counts or future scenario estimates. The original inputs were not reprocessed here. Independent station validation is still needed to assess local temperature accuracy; coarse grids can smooth coastal and mountain extremes. Findings apply to the study rectangle, and mean-duration estimates remain conditional on event-bearing years.

## Run

Open a terminal in this folder, with R installed:

```sh
Rscript install_dependencies.R
Rscript run_analysis.R
```

In RStudio, set the working directory to this folder and run:

```r
source("install_dependencies.R")
source("run_analysis.R")
```

The default analysis reads the **bundled, compressed NASA POWER responses** and runs without downloading the climate data again. A local RDS cache speeds later runs. Changing the area or analysis years selects a different cache/snapshot key and can require API access. Change the settings at the top of `Heat.R`; keep the full reference period inside the analysis period.

This release was checked with R 4.3.3, jsonlite 1.8.8 and ggplot2 3.4.4; the actual package session is recorded in `results/sessionInfo.txt`.

`Heat.R` runs the primary analysis and bandwidth/common-year checks. `scripts/check_minimum_years.R` reuses its functions and annual indices to evaluate 10, 15 and 20 available-year rules. `run_analysis.R` runs both.

## Data and methods

| Choice | Specification |
| --- | --- |
| Input | NASA POWER `T2M_MAX`, degrees Celsius, UTC; 1981–2025 |
| Domain | 36–42°E, 32.5–37°N; 100 native cells, 0.625° longitude × 0.5° latitude |
| Reference | Fixed 1991–2020; calendar-day P90 and P95 in circular ±2-day windows; quantile type 8 |
| Leap day | Excluded from quantile fitting; threshold interpolated from 28 February and 1 March; observed 29 February retained in events |
| Event | A maximal consecutive run with Tmax **strictly above** its threshold, lasting at least 2, 3 or 4 days |
| Annual indices | Number of events, days in qualifying events, and conditional mean event duration |
| Primary change | OLS slope per decade; Sen slope as a descriptive robustness estimate |
| Uncertainty | Approximate pointwise HAC intervals and tests; Bartlett weights use actual calendar-year gaps |
| Eligibility | At least 15 available annual values; a pragmatic sparse-sample guard, checked against 10 and 20 |
| Multiplicity | BY at FDR 0.05 across all 1,800 planned primary cases, including unestimable cases |

Each cell must have exactly the expected daily dates, no duplicates or missing temperatures, and the expected units/time standard. Incomplete input stops the analysis instead of being silently imputed.

Events are detected before grouping by year. Event counts and full event durations belong to the **onset year**; event days belong to their **actual calendar year**. Years affected by a record-edge run are conservatively excluded. A genuine no-event year has count/day totals of zero and mean duration `NA`.

Conditional duration uses each definition's available event-bearing years. Common event-years are a secondary comparison. The primary HAC lag is 3 years, using `floor(4 * (45/100)^(2/9))`; lags 2 and 6 are retained as checks. Covariance uses the finite-sample factor `n/(n-2)` and t critical values with `n-2` degrees of freedom. These are approximate intervals, conditional on the fixed thresholds and product.

The rectangle is the original extraction domain, **not a Syria boundary or a land-only mask**. Results apply to the returned cells. Displayed cells are clipped to the requested rectangle, so the northern and southern edge rows appear at half height; their calculated values are unchanged.

## Outputs

All outputs are written to `results/`:

- `heatwave_annual.csv`: one row per cell, year and definition.
- `heatwave_trends.csv`: primary slopes, coverage, pointwise intervals, p-values and BY adjustments.
- `heatwave_trend_checks.csv`: alternative HAC lags and common-year duration fits; exploratory checks without an additional adjusted testing family.
- `minimum_years_trends.csv` and `minimum_years_sensitivity.csv`: complete eligibility-rule comparisons and their summary.
- Five PNG figures: three trend maps, duration coverage, and minimum-year sensitivity.
- `analysis_settings.rds` and `sessionInfo.txt`: settings and the actual R session.

Slope units are **events/year per decade** for frequency, **event-days/year per decade** for day totals, and **days/event per decade** for conditional mean duration. Grid cells and nested definitions are correlated; counts of positive or significant cases are descriptive summaries.

## Scope and checks

The final workflow was run on the bundled full-period snapshot. Counts and mean durations were matched to the earlier implementation for all 27,000 annual rows. Calendar-year event days were independently checked against event/year interval overlaps. HAC calculations were checked against `sandwich::NeweyWest` for complete samples and an independent calendar-gap calculation for duration samples. The 10/15/20-year comparisons are reproducible in R.

POWER is a gridded reanalysis-derived product and can smooth local extremes, particularly in complex coastal/topographic settings. These full-year relative spells include winter warm anomalies; they are not health-risk thresholds or standard six-day WSDI. The intervals exclude threshold-estimation and product uncertainty, and sparse duration estimates remain approximate. No station-based validation, climate attribution or future scenario modelling is performed here.

## Data provenance, references and licence

The 45 annual API responses are archived in `data/power_snapshot/` as gzip JSON. `data/snapshot_manifest.csv` records SHA-256 hashes of both the archived files and their original JSON, API version and time standard. The archive was packaged for this release on 2026-10-25; the original download timestamp was not retained. Cache files are excluded from Git.

Data acknowledgement: **NASA POWER Project**, supported through NASA's Applied Sciences Program at NASA Langley Research Center. The included responses identify MERRA2 and POWER as sources and POWER Daily API v2.10.1. The MIT licence in `LICENSE` applies to project code; data attribution and NASA's source terms remain applicable.

- [NASA POWER meteorological data and native resolution](https://power.larc.nasa.gov/docs/methodology/meteorology/)
- [NASA POWER daily API, including UTC requests](https://power.larc.nasa.gov/docs/services/api/temporal/daily/)
- [NOAA GHCNd station-data overview](https://www.ncei.noaa.gov/products/land-based-station/global-historical-climatology-network-daily)
- [Berkeley Earth station-based reconstruction methodology](https://berkeleyearth.org/methodology/)
- [Climdex percentile/spell index definitions](https://www.climdex.org/learn/indices/)
- [Newey–West estimator documentation](https://search.r-project.org/CRAN/refmans/sandwich/html/NeweyWest.html)
- [R multiple-testing adjustments](https://stat.ethz.ch/R-manual/R-patched/library/stats/html/p.adjust.html)
