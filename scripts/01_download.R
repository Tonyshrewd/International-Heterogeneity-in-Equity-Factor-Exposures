# =============================================================================
# 01_download.R  --  Download JKP Global Factor Data (country level, monthly)
# FIN 864 project: International Heterogeneity in Equity Factor Exposures
#
# Run from the project root:   Rscript scripts/01_download.R
#
# Downloads one file per dataset covering ALL countries:
#   mkt          country market returns                (dependent variable)
#   all_themes   13 theme (cluster) returns             (main regressors)
#   all_factors  153 individual factor returns          (within-theme tests)
# in three weighting schemes:
#   vw_cap  capped value weighted  (main specification)
#   vw      value weighted         (robustness)
#   ew      equal weighted         (robustness)
#
# The server's Last-Modified date is saved as the data release stamp in
# raw/download_log.csv, as required for the feasibility report.
# Re-running skips files that are already fully downloaded.
# =============================================================================

base_url   <- "https://jkpfactors-data.s3.amazonaws.com/public/"
location   <- "all_countries"
datasets   <- c("mkt", "all_themes", "all_factors")
weightings <- c("vw_cap", "vw", "ew")
frequency  <- "monthly"

dir.create("raw/zips", recursive = TRUE, showWarnings = FALSE)
dir.create("raw/csv",  recursive = TRUE, showWarnings = FALSE)
options(timeout = 1800)  # large files: allow up to 30 minutes each

# File names on the server look like  [all_countries]_[mkt]_[monthly]_[vw_cap].zip
make_url <- function(loc, ds, freq, w) {
  file <- sprintf("[%s]_[%s]_[%s]_[%s].zip", loc, ds, freq, w)
  paste0(base_url, URLencode(file, reserved = TRUE))
}

get_header <- function(headers, field) {
  line <- grep(paste0("^", field, ":"), headers, ignore.case = TRUE, value = TRUE)
  if (length(line) == 0) return(NA_character_)
  trimws(sub("^[^:]+:", "", line[length(line)]))
}

log_rows <- list()

for (w in weightings) {
  for (ds in datasets) {
    url  <- make_url(location, ds, frequency, w)
    dest <- file.path("raw/zips", sprintf("%s_%s_%s_%s.zip", location, ds, frequency, w))

    hdr      <- curlGetHeaders(url)
    status   <- attr(hdr, "status")
    if (status != 200) stop("Server returned status ", status, " for ", url)
    size     <- as.numeric(get_header(hdr, "Content-Length"))
    last_mod <- get_header(hdr, "Last-Modified")

    if (file.exists(dest) && file.size(dest) == size) {
      message("Already downloaded: ", basename(dest))
    } else {
      message("Downloading ", basename(dest), " (", round(size / 1e6, 1), " MB) ...")
      download.file(url, dest, mode = "wb", quiet = TRUE)
      if (file.size(dest) != size) stop("Incomplete download: ", dest)
    }

    csv_name <- unzip(dest, list = TRUE)$Name
    unzip(dest, exdir = "raw/csv", overwrite = TRUE)

    log_rows[[length(log_rows) + 1]] <- data.frame(
      dataset       = ds,
      weighting     = w,
      frequency     = frequency,
      location      = location,
      url           = url,
      zip_file      = dest,
      csv_file      = file.path("raw/csv", csv_name),
      bytes         = size,
      server_last_modified = last_mod,
      downloaded_on = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    )
  }
}

log <- do.call(rbind, log_rows)
write.csv(log, "raw/download_log.csv", row.names = FALSE)
message("\nDone. ", nrow(log), " files. Release stamp(s): ",
        paste(unique(log$server_last_modified), collapse = "; "))
