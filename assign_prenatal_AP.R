################################################################################
# Assign prenatal air quality exposures (PM2.5, O3, NO2) to subjects
#
# Data source:
#   NASA/SEDAC "Daily and Annual PM2.5, O3 and NO2 Concentrations at ZIP Codes
#   for the Contiguous U.S., 2000-2016"
#   https://data.nasa.gov/dataset/daily-and-annual-pm2-5-o3-and-no2-concentrations-at-zip-codes-for-the-contiguous-u-s-2000-
#
# Method:
#   For each subject, compute the mean daily concentration of PM2.5, O3, and
#   NO2 over the 9 calendar months immediately preceding their date of birth
#   (i.e., from the first day of [birth_month - 9] to DOB - 1 day, inclusive).
#
# Expected input files:
#   1. subjects.csv  — subject-level data with at minimum:
#        subject_id  : unique subject identifier
#        dob         : date of birth (any format parseable by as.Date / lubridate)
#        zip         : 5-digit ZIP code of residence during pregnancy
#
#   2. Air-quality CSV(s) downloaded from the NASA data portal above.
#      Expected columns (rename below if yours differ):
#        ZIP         : 5-digit ZIP code
#        Date        : date in YYYY-MM-DD format  (a.k.a. "date" in some versions)
#        PM25        : daily PM2.5 concentration (µg/m³)
#        O3          : daily O3 concentration (ppb)
#        NO2         : daily NO2 concentration (ppb)
#
#      The script accepts a single CSV or a directory of CSVs (one per year is
#      the common distribution format).  Set AQ_INPUT below accordingly.
################################################################################

library(dplyr)
library(lubridate)
library(readr)
library(purrr)

# ── USER CONFIGURATION ────────────────────────────────────────────────────────

SUBJECTS_FILE <- "subjects.csv"   # path to subject-level file

# Path to air-quality data: either a single CSV or a directory of CSVs
AQ_INPUT      <- "aq_data/"       # e.g. "aq_data/" or "aq_data/daily_2010.csv"

OUTPUT_FILE   <- "subjects_with_aq_exposure.csv"

# Column name mappings — edit if your downloaded files use different names
AQ_COLS <- list(
  zip  = "ZIP",
  date = "Date",
  pm25 = "PM25",
  o3   = "O3",
  no2  = "NO2"
)

# ── LOAD DATA ─────────────────────────────────────────────────────────────────

message("Loading subject data...")
subjects <- read_csv(SUBJECTS_FILE, show_col_types = FALSE) |>
  mutate(dob = as.Date(dob),
         zip = as.character(zip) |> stringr::str_pad(5, pad = "0"))

message("Loading air-quality data...")
load_aq <- function(path) {
  if (dir.exists(path)) {
    csv_files <- list.files(path, pattern = "\\.csv$", full.names = TRUE,
                            recursive = TRUE)
    if (length(csv_files) == 0) stop("No CSV files found in directory: ", path)
    map(csv_files, \(f) read_csv(f, show_col_types = FALSE)) |> list_rbind()
  } else {
    read_csv(path, show_col_types = FALSE)
  }
}

aq_raw <- load_aq(AQ_INPUT)

# Standardise column names
aq <- aq_raw |>
  rename(
    zip  = all_of(AQ_COLS$zip),
    date = all_of(AQ_COLS$date),
    pm25 = all_of(AQ_COLS$pm25),
    o3   = all_of(AQ_COLS$o3),
    no2  = all_of(AQ_COLS$no2)
  ) |>
  mutate(
    zip  = as.character(zip) |> stringr::str_pad(5, pad = "0"),
    date = as.Date(date)
  ) |>
  select(zip, date, pm25, o3, no2) |>
  filter(!is.na(date))

message(sprintf("Air-quality data: %d rows spanning %s to %s",
                nrow(aq), min(aq$date), max(aq$date)))

# ── COMPUTE PRENATAL EXPOSURE WINDOWS ─────────────────────────────────────────
# Prenatal window: 9 calendar months before DOB
#   window_end   = DOB - 1 day
#   window_start = first day of (DOB month - 9 months)
#                = floor_date(DOB %m-% months(9), "month")

subjects <- subjects |>
  mutate(
    window_end   = dob - days(1),
    window_start = floor_date(dob %m-% months(9), unit = "month")
  )

message("Prenatal window example (first subject):")
message(sprintf("  DOB: %s  |  window: %s to %s",
                subjects$dob[1], subjects$window_start[1],
                subjects$window_end[1]))

# ── ASSIGN EXPOSURE VALUES ─────────────────────────────────────────────────────
# For each subject, filter the AQ data to their ZIP + date window, then average.
# Uses a join + filter strategy that scales well for large cohorts.

message("Assigning exposure values...")

# Create a lookup key: subject x (zip, window_start, window_end)
subject_windows <- subjects |>
  select(subject_id, zip, window_start, window_end)

# Join AQ to subject windows on ZIP, then filter to date window
exposure <- subject_windows |>
  left_join(aq, by = "zip", relationship = "many-to-many") |>
  filter(date >= window_start & date <= window_end) |>
  group_by(subject_id) |>
  summarise(
    pm25_mean   = mean(pm25, na.rm = TRUE),
    o3_mean     = mean(o3,   na.rm = TRUE),
    no2_mean    = mean(no2,  na.rm = TRUE),
    n_days_pm25 = sum(!is.na(pm25)),
    n_days_o3   = sum(!is.na(o3)),
    n_days_no2  = sum(!is.na(no2)),
    .groups = "drop"
  )

# ── MERGE BACK AND REPORT ──────────────────────────────────────────────────────

result <- subjects |>
  left_join(exposure, by = "subject_id")

n_missing <- sum(is.na(result$pm25_mean))
if (n_missing > 0) {
  warning(sprintf(
    "%d subject(s) have no matching AQ data. Check ZIP codes and date ranges.",
    n_missing))
}

message(sprintf("\nExposure summary (n = %d subjects):", nrow(result)))
message(sprintf("  PM2.5: mean=%.2f, range=[%.2f, %.2f] µg/m³",
                mean(result$pm25_mean, na.rm = TRUE),
                min(result$pm25_mean,  na.rm = TRUE),
                max(result$pm25_mean,  na.rm = TRUE)))
message(sprintf("  O3:    mean=%.2f, range=[%.2f, %.2f] ppb",
                mean(result$o3_mean,   na.rm = TRUE),
                min(result$o3_mean,    na.rm = TRUE),
                max(result$o3_mean,    na.rm = TRUE)))
message(sprintf("  NO2:   mean=%.2f, range=[%.2f, %.2f] ppb",
                mean(result$no2_mean,  na.rm = TRUE),
                min(result$no2_mean,   na.rm = TRUE),
                max(result$no2_mean,   na.rm = TRUE)))

write_csv(result, OUTPUT_FILE)
message(sprintf("\nResults written to: %s", OUTPUT_FILE))
