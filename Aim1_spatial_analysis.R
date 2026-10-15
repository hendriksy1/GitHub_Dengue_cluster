# Aim 1: space-time scan of laboratory-confirmed dengue cases, Valle del Cauca, 2022
# R packages: terra, readxl, lubridate. SaTScan 10.3.2

# File paths
case_file       <- "path/to/surveillance_cases.xlsx"
population_file <- "path/to/population_raster.tif"
satscan_binary  <- "path/to/satscan"
output_dir      <- "path/to/private_outputs/spatial"

municipalities <- c("Cali", "Palmira", "Yumbo", "Jamundi", "Buga", "Tulua", "Florida", "Cartago", "Candelaria")
maximum_months <- c(6L, 3L, 1L)

# Rectangles for case selection: c(lat_min, lat_max, lon_min, lon_max)
bounds <- list(
  Candelaria = c(min(3.396209, 3.400042, 3.423978, 3.425073), max(3.396209, 3.400042, 3.423978, 3.425073), min(-76.360572, -76.356497, -76.342941, -76.341452), max(-76.360572, -76.356497, -76.342941, -76.341452)),
  Cali = c(min(3.285969, 3.302464, 3.513461, 3.514030), max(3.285969, 3.302464, 3.513461, 3.514030), min(-76.454641, -76.586819, -76.447804, -76.595365), max(-76.454641, -76.586819, -76.447804, -76.595365)),
  Palmira = c(min(3.498942, 3.499456, 3.562165, 3.565420), max(3.498942, 3.499456, 3.562165, 3.565420), min(-76.263677, -76.324273, -76.263505, -76.325990), max(-76.263677, -76.324273, -76.263505, -76.325990)),
  Yumbo = c(3.542880, 3.604621, -76.517674, -76.474992),
  Jamundi = c(min(3.241715, 3.240693, 3.280452, 3.281986), max(3.241715, 3.240693, 3.280452, 3.281986), min(-76.560897, -76.558848, -76.529012, -76.525299), max(-76.560897, -76.558848, -76.529012, -76.525299)),
  Buga = c(min(3.857909, 3.861028, 3.937103, 3.936715), max(3.857909, 3.861028, 3.937103, 3.936715), min(-76.284782, -76.322462, -76.276306, -76.316469), max(-76.284782, -76.322462, -76.276306, -76.316469)),
  Tulua = c(min(4.060501, 4.060329, 4.122312, 4.123853), max(4.060501, 4.060329, 4.122312, 4.123853), min(-76.171906, -76.224091, -76.171048, -76.228555), max(-76.171906, -76.224091, -76.171048, -76.228555)),
  Florida = c(min(3.305268, 3.306291, 3.341061, 3.343106), max(3.305268, 3.306291, 3.341061, 3.343106), min(-76.257840, -76.253059, -76.208556, -76.207304), max(-76.257840, -76.253059, -76.208556, -76.207304)),
  Cartago = c(min(4.718791, 4.724292, 4.772476, 4.771907), max(4.718791, 4.724292, 4.772476, 4.771907), min(-75.897838, -75.957989, -75.893079, -75.961606), max(-75.897838, -75.957989, -75.893079, -75.961606))
)

library(terra)
library(readxl)
library(lubridate)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
input_dir <- file.path(output_dir, "inputs")
run_name <- function(m) paste0("Satscan_", tolower(m), "_20_10cases_9999")
input_file <- function(m, ext) file.path(input_dir, m, paste0(run_name(m), ext))
result_stem <- function(m, maximum) file.path(output_dir, paste0("results_", maximum, "month"), m, run_name(m))
write_ascii <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ", na = "")
}
numeric_value <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
parse_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  raw <- trimws(as.character(x))
  raw[raw %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  result <- rep(as.Date(NA), length(raw))
  number <- suppressWarnings(as.numeric(raw))
  serial <- !is.na(number) & number >= 20000 & number <= 60000
  result[serial] <- as.Date(number[serial], origin = "1899-12-30")
  rest <- is.na(result) & !is.na(raw)
  result[rest] <- as.Date(lubridate::parse_date_time(raw[rest],
                                                     orders = c("ymd", "dmy", "mdy", "ymd HMS", "dmy HMS", "mdy HMS"), tz = "UTC", quiet = TRUE))
  result
}

# 1. Laboratory confirmation and consultation dates
all_cases <- as.data.frame(readxl::read_excel(case_file))
stopifnot(all(c("fec_con_", "tip_cas_", "ajuste_", "latitude", "longitude") %in% names(all_cases)))
all_cases$date <- parse_date(all_cases$fec_con_)
all_cases$latitude <- numeric_value(all_cases$latitude)
all_cases$longitude <- numeric_value(all_cases$longitude)
lab <- all_cases[which(all_cases$tip_cas_ == 3 | all_cases$ajuste_ == 3), ]
print(table(tip_cas = lab$tip_cas_, ajuste = lab$ajuste_, useNA = "ifany"))
population_full <- terra::rast(population_file)
stopifnot(terra::is.lonlat(population_full), terra::nlyr(population_full) == 1L)

# 2. SaTScan input files by municipality
flow <- data.frame(municipality = municipalities, in_extent = NA_integer_, analyzed = NA_integer_)
for (m in municipalities) {
  b <- bounds[[m]]
  d <- lab[which(lab$latitude >= b[1] & lab$latitude <= b[2] &
                   lab$longitude >= b[3] & lab$longitude <= b[4]), ]
  flow$in_extent[flow$municipality == m] <- nrow(d)
  # No scan where the selection contains no cases
  if (nrow(d) == 0L) {
    flow$analyzed[flow$municipality == m] <- 0L
    next
  }
  stopifnot(!anyNA(d$date), all(format(d$date, "%Y") == "2022"))
  # Population raster cropped to the case-point extent
  points <- terra::vect(d, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  points <- terra::project(points, terra::crs(population_full))
  grid <- terra::crop(population_full, points)
  xy <- terra::xyFromCell(grid, seq_len(terra::ncell(grid)))
  pop <- data.frame(ID = seq_len(nrow(xy)), pop = as.vector(terra::values(grid)))
  geo <- data.frame(ID = pop$ID, Latitude = xy[, 2], Longitude = xy[, 1])
  populated <- !is.na(pop$pop) & pop$pop > 0
  pop <- pop[populated, ]
  geo <- geo[populated, ]
  cell <- terra::cellFromXY(grid, terra::crds(points))
  keep <- !is.na(cell)
  stopifnot(all(cell[keep] %in% pop$ID))
  flow$analyzed[flow$municipality == m] <- sum(keep)
  d <- d[keep, ]
  cell <- cell[keep]
  counts <- aggregate(list(count = rep.int(1L, nrow(d))), list(ID = cell, date = d$date), sum)
  counts <- counts[order(counts$ID, counts$date), ]
  write_ascii(data.frame(ID = counts$ID, cases = counts$count, date = format(counts$date, "%Y/%m/%d")), input_file(m, ".cas"))
  write_ascii(geo, input_file(m, ".geo"))
  write_ascii(data.frame(ID = pop$ID, census_time = "2022/1/1", population = pop$pop), input_file(m, ".pop"))
}
scanned_municipalities <- flow$municipality[flow$analyzed > 0L]
case_flow <- data.frame(
  step = c("notified cases", "laboratory-confirmed", "analyzed in spatial scans"),
  n = c(nrow(all_cases), nrow(lab), sum(flow$analyzed)))
print(case_flow)
print(flow)
write.csv(case_flow, file.path(output_dir, "case_flow.csv"), row.names = FALSE)
write.csv(flow, file.path(output_dir, "cases_by_municipality.csv"), row.names = FALSE)

# 3. Poisson space-time scans
# Spatial limits: 20% of population, 1 km radius, at least 10 cases
# Monthly aggregation, 9,999 replications, maximum windows of 6, 3 and 1 months
parameters <- c(
  "[Input]", "CaseFile=", "PopulationFile=", "CoordinatesFile=",
  "PrecisionCaseTimes=2", "StartDate=2022/1/1", "EndDate=2022/12/31", "CoordinatesType=1",
  "[Analysis]", "AnalysisType=3", "ModelType=0", "ScanAreas=1",
  "TimeAggregationUnits=2", "TimeAggregationLength=1",
  "[Spatial Window]", "MaxSpatialSizeInPopulationAtRisk=20",
  "UseDistanceFromCenterOption=y", "MaxSpatialSizeInDistanceFromCenter=1",
  "IncludePurelyTemporal=n", "SpatialWindowShapeType=0",
  "[Temporal Window]", "MinimumTemporalClusterSize=1",
  "MaxTemporalSizeInterpretation=0", "MaxTemporalSize=50", "IncludePurelySpatial=n",
  "[Cluster Restrictions]", "MinimumCasesInHighRateClusters=10",
  "[Space and Time Adjustments]", "TimeTrendAdjustmentType=0", "SpatialAdjustmentType=0",
  "[Inference]", "PValueReportType=0", "MonteCarloReps=9999", "IterativeScan=n",
  "[Output]", "ResultsFile=", "MostLikelyClusterEachCentroidASCII=y",
  "CensusAreasReportedClustersASCII=y", "OutputGoogleEarthKML=n", "OutputGoogleMaps=n",
  "OutputShapefiles=n", "OutputCartesianGraph=n",
  "[Spatial Output]", "CriteriaForReportingSecondaryClusters=0",
  "UseReportOnlySmallerClusters=y", "MaxSpatialSizeInPopulationAtRisk_Reported=20",
  "UseDistanceFromCenterOption_Reported=y", "MaxSpatialSizeInDistanceFromCenter_Reported=1",
  "[Other Output]", "PrintAsciiColumnHeaders=n",
  "[Cluster Drilldown]", "PerformStandardDrilldown=n",
  "[Run Options]", "NumberParallelProcesses=4", "LogRunToHistoryFile=n",
  "[System]", "Version=10.3.2", "RandomlyGenerateSeed=n", "RandomSeed=20260622"
)
set_parameter <- function(x, key, value) {
  i <- grep(paste0("^", key, "="), x)
  x[i] <- paste0(key, "=", value)
  x
}
for (maximum in maximum_months) for (m in scanned_municipalities) {
  stem <- result_stem(m, maximum)
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  stem <- file.path(normalizePath(dirname(stem)), basename(stem))
  prm <- parameters
  prm <- set_parameter(prm, "CaseFile", normalizePath(input_file(m, ".cas")))
  prm <- set_parameter(prm, "PopulationFile", normalizePath(input_file(m, ".pop")))
  prm <- set_parameter(prm, "CoordinatesFile", normalizePath(input_file(m, ".geo")))
  prm <- set_parameter(prm, "ResultsFile", stem)
  if (maximum != 6L) {
    prm <- set_parameter(prm, "MaxTemporalSizeInterpretation", "1")
    prm <- set_parameter(prm, "MaxTemporalSize", as.character(maximum))
  }
  writeLines(prm, paste0(stem, ".prm"))
  log_file <- paste0(stem, ".run.log")
  message(m, ": maximum ", maximum, " months")
  status <- system2(satscan_binary, shQuote(paste0(stem, ".prm")), stdout = log_file, stderr = log_file)
  if (status != 0L || !any(grepl("SaTScan completed successfully.", readLines(log_file, warn = FALSE), fixed = TRUE)))
    stop("SaTScan did not complete: ", log_file)
}

# 4. Clusters with p < 0.05 and their location-month membership
read_ascii <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (!any(nzchar(trimws(lines)))) return(data.frame())
  read.table(text = lines, header = FALSE, stringsAsFactors = FALSE, colClasses = "character")
}
case_counts <- function(m) {
  x <- read_ascii(input_file(m, ".cas"))
  dates <- as.Date(x[[3]], "%Y/%m/%d")
  tapply(as.numeric(x[[2]]), paste(x[[1]], as.integer(format(dates, "%m")), sep = ":"), sum)
}
count_cells <- function(counts, cells) sum(counts[intersect(names(counts), cells)])
membership <- function(clusters, key) unique(unlist(lapply(clusters, `[[`, key), use.names = FALSE))
read_clusters <- function(m, maximum, counts) {
  stem <- result_stem(m, maximum)
  x <- read_ascii(paste0(stem, ".col.txt"))
  gis <- read_ascii(paste0(stem, ".gis.txt"))
  if (!nrow(x)) return(list())
  columns <- c("cluster", "center_location_id", "latitude", "longitude", "radius_km", "span_km",
               "start_date", "end_date", "n_locations", "log_likelihood_ratio", "p_value",
               "observed_cases", "expected_cases", "observed_expected", "relative_risk", "population", "flag")
  if (ncol(x) == 16L) columns <- setdiff(columns, "span_km")
  stopifnot(ncol(x) == length(columns), ncol(gis) == 16L)
  names(x) <- columns
  for (key in setdiff(columns, c("start_date", "end_date", "flag"))) x[[key]] <- as.numeric(x[[key]])
  x <- x[x$p_value < 0.05, , drop = FALSE]
  lapply(seq_len(nrow(x)), function(i) {
    d <- as.list(x[i, , drop = FALSE])
    d$municipality <- m
    d$maximum_months <- maximum
    first <- as.Date(d$start_date, "%Y/%m/%d")
    last <- as.Date(d$end_date, "%Y/%m/%d")
    months <- seq(as.integer(format(first, "%m")), as.integer(format(last, "%m")))
    locations <- unique(gis[[1]][as.integer(gis[[2]]) == d$cluster])
    cells <- as.vector(outer(locations, months, paste, sep = ":"))
    stopifnot(length(locations) == d$n_locations, length(months) <= maximum,
              count_cells(counts, cells) == d$observed_cases)
    d$months <- months
    d$cells <- cells
    d
  })
}
counts <- setNames(lapply(scanned_municipalities, case_counts), scanned_municipalities)
clusters <- list()
for (m in scanned_municipalities) for (maximum in maximum_months)
  clusters[[paste(m, maximum)]] <- read_clusters(m, maximum, counts[[m]])

# 5. Primary results and temporal sensitivity analysis
# No eligible cases: scan status "Not scanned", cluster count NA
season <- function(months) {
  dry <- any(months %in% c(12, 1, 2, 6, 7, 8, 9))
  wet <- any(months %in% c(3, 4, 5, 10, 11))
  if (dry && wet) "Both" else if (dry) "Warm-dry" else "Cool-wet"
}
cluster_results <- do.call(rbind, lapply(names(clusters), function(k) {
  lapply(clusters[[k]], function(c) data.frame(
    municipality = c$municipality, maximum_months = c$maximum_months, cluster = c$cluster,
    latitude = c$latitude, longitude = c$longitude, radius_km = c$radius_km,
    first_month = min(c$months), last_month = max(c$months), duration_months = length(c$months),
    observed_cases = c$observed_cases, population = c$population,
    relative_risk = c$relative_risk, p_value = c$p_value, season = season(c$months),
    reached_maximum = length(c$months) == c$maximum_months))
}) |> unlist(recursive = FALSE))
write.csv(cluster_results, file.path(output_dir, "cluster_results.csv"), row.names = FALSE)

iqr <- function(x) quantile(x, c(.25, .5, .75), names = FALSE)
table1 <- do.call(rbind, lapply(municipalities, function(m) {
  scanned <- m %in% scanned_municipalities
  z <- clusters[[paste(m, 6L)]]
  cr <- cluster_results[cluster_results$municipality == m & cluster_results$maximum_months == 6L, ]
  inside <- if (scanned) count_cells(counts[[m]], membership(z, "cells")) else NA_real_
  n <- flow$analyzed[flow$municipality == m]
  q <- function(v) if (nrow(cr)) iqr(cr[[v]]) else rep(NA_real_, 3)
  data.frame(municipality = m, first_month = if (nrow(cr)) min(cr$first_month) else NA,
             last_month = if (nrow(cr)) max(cr$last_month) else NA,
             scan_status = if (scanned) "Scanned" else "Not scanned", clusters = if (scanned) nrow(cr) else NA_integer_,
             cases_in_study_area = n, cases_in_clusters = inside, percent_in_clusters = round(100 * inside / n, 1),
             population_q1 = q("population")[1], population_median = q("population")[2], population_q3 = q("population")[3],
             radius_q1 = q("radius_km")[1], radius_median = q("radius_km")[2], radius_q3 = q("radius_km")[3],
             cases_q1 = q("observed_cases")[1], cases_median = q("observed_cases")[2], cases_q3 = q("observed_cases")[3],
             rr_q1 = q("relative_risk")[1], rr_median = q("relative_risk")[2], rr_q3 = q("relative_risk")[3])
}))
write.csv(table1, file.path(output_dir, "table1_primary_scan.csv"), row.names = FALSE)

# Median RR rounded to two decimals, with halves rounded up
s3 <- do.call(rbind, lapply(maximum_months, function(maximum) {
  cr <- cluster_results[cluster_results$maximum_months == maximum, ]
  unique_n <- sum(sapply(scanned_municipalities, function(m) count_cells(counts[[m]], membership(clusters[[paste(m, maximum)]], "cells"))))
  retained <- sum(sapply(scanned_municipalities, function(m) count_cells(counts[[m]],
                                                                         intersect(membership(clusters[[paste(m, maximum)]], "cells"), membership(clusters[[paste(m, 6L)]], "cells")))))
  primary_n <- sum(sapply(scanned_municipalities, function(m) count_cells(counts[[m]], membership(clusters[[paste(m, 6L)]], "cells"))))
  total <- sum(sapply(counts, sum))
  data.frame(maximum_months = maximum, cases_analyzed = total, clusters = nrow(cr),
             unique_cases = unique_n, percent_cases = round(100 * unique_n / total, 1),
             primary_cases_retained = retained, percent_primary_retained = round(100 * retained / primary_n, 1),
             duration_median = median(cr$duration_months), duration_min = min(cr$duration_months), duration_max = max(cr$duration_months),
             reached_maximum = sum(cr$reached_maximum), percent_reached_maximum = round(100 * mean(cr$reached_maximum), 1),
             rr_median = floor(median(cr$relative_risk) * 100 + 0.5) / 100, radius_min = min(cr$radius_km), radius_max = max(cr$radius_km),
             warm_dry_only = sum(cr$season == "Warm-dry"), cool_wet_only = sum(cr$season == "Cool-wet"), both_seasons = sum(cr$season == "Both"),
             t(sapply(municipalities, function(m) if (m %in% scanned_municipalities) sum(cr$municipality == m) else NA_integer_) |> setNames(paste0("n_", municipalities))))
}))
print(table1, row.names = FALSE)
print(s3, row.names = FALSE)
write.csv(s3, file.path(output_dir, "s3_table_temporal_windows.csv"), row.names = FALSE)

# Software versions used for the analysis.
sessionInfo()
