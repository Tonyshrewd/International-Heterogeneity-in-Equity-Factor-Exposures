# =============================================================================
# 02_build_panel.R  --  Clean JKP data and build the country-month panel
# FIN 864 project: International Heterogeneity in Equity Factor Exposures
#
# Run from the project root (after 01_download.R):
#   Rscript scripts/02_build_panel.R
#
# Inputs
#   raw/csv/[all_countries]_[<dataset>]_[monthly]_[<weighting>].csv
#   raw/factor_details.xlsx        official JKP factor list (names, citations)
#   scripts/theme_map.csv          factor -> theme, from Table 9 of the JKP
#                                  documentation (the xlsx has no theme column)
#   scripts/country_meta.csv       country names and MSCI classification
#                                  (Table 10 of the JKP documentation)
#
# Outputs (data/)
#   panel_<w>.csv / .rds           one row per country-month: market return,
#                                  13 theme returns  (w = vw_cap, vw, ew)
#   market_long.csv                market returns, all weightings, long format
#   themes_long.csv                theme returns, all weightings, long format
#   factors_long_<w>.csv.gz/.rds   153 factor returns with theme labels
#   factor_dictionary.csv          factor names, themes, signs, citations
# Outputs (output/)
#   check_theme_reconstruction.csv verifies theme = average of its factors
#   coverage_by_country.csv        sample construction (Table 1 input)
#   coverage_theme_by_country.csv  months available per theme and country
#   common_window_feasibility.csv  common-history window, 5 pilot countries
#   checks_log.txt                 all data checks, in plain text
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

weightings  <- c("vw_cap", "vw", "ew")
w_main      <- "vw_cap"
pilot       <- c("USA", "JPN", "DEU", "GBR", "AUS")
min_months  <- 120

dir.create("data",   showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)

log_file <- "output/checks_log.txt"
cat("JKP data build log --", format(Sys.time()), "\n\n", file = log_file)
say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n")
  cat(txt, "\n", file = log_file, append = TRUE)
}

read_jkp <- function(ds, w) {
  f <- sprintf("raw/csv/[all_countries]_[%s]_[monthly]_[%s].csv", ds, w)
  if (!file.exists(f)) stop("Missing file: ", f, "  -- run 01_download.R first")
  dt <- fread(f, showProgress = FALSE)
  dt[, `:=`(country = toupper(location), date = as.IDate(date))]
  dt[]
}

# -----------------------------------------------------------------------------
# 1. Reference tables
# -----------------------------------------------------------------------------
theme_map <- fread("scripts/theme_map.csv")
themes    <- sort(unique(theme_map$theme))
countries <- fread("scripts/country_meta.csv", encoding = "UTF-8")

details <- as.data.table(read_excel("raw/factor_details.xlsx"))[!is.na(abr_jkp)]
factor_dictionary <- merge(
  theme_map,
  details[, .(factor = abr_jkp, description = name_new, citation = cite,
              original_sample = `in-sample period`, direction_xlsx = direction,
              significant_in_original = significance, hxz_group = group)],
  by = "factor", all.x = TRUE
)
stopifnot(nrow(factor_dictionary) == 153, !anyNA(factor_dictionary$description))
say("Reference tables: ", nrow(theme_map), " factors in ", length(themes), " themes")

# -----------------------------------------------------------------------------
# 2. Read and stack the three datasets for each weighting
# -----------------------------------------------------------------------------
mkt_list <- th_list <- list()
fac_main <- NULL

for (w in weightings) {
  m <- read_jkp("mkt", w)
  t <- read_jkp("all_themes", w)

  # Key checks: one observation per country-month (and theme)
  stopifnot(!anyDuplicated(m, by = c("country", "date")))
  stopifnot(!anyDuplicated(t, by = c("country", "date", "name")))
  stopifnot(setequal(unique(t$name), themes))

  mkt_list[[w]] <- m[, .(country, date, weighting = w, mkt = ret, mkt_n_stocks = n_stocks)]
  th_list[[w]]  <- t[, .(country, date, weighting = w, theme = name, ret, n_factors)]

  say(sprintf("[%s] market: %d rows, %d countries | themes: %d rows, %d countries",
              w, nrow(m), uniqueN(m$country), nrow(t), uniqueN(t$country)))

  # Factor file: keep all weightings on disk, but only the main one in memory
  f <- read_jkp("all_factors", w)
  stopifnot(!anyDuplicated(f, by = c("country", "date", "name")))
  stopifnot(setequal(unique(f$name), theme_map$factor))
  f <- merge(f[, .(country, date, factor = name, direction, n_stocks, n_stocks_min, ret)],
             theme_map, by = "factor", all.x = TRUE)
  setcolorder(f, c("country", "date", "theme", "factor"))
  setorder(f, country, theme, factor, date)
  fwrite(f, sprintf("data/factors_long_%s.csv.gz", w))
  saveRDS(f, sprintf("data/factors_long_%s.rds", w))
  say(sprintf("[%s] factors: %d rows, %d countries, %d factors",
              w, nrow(f), uniqueN(f$country), uniqueN(f$factor)))
  if (w == w_main) fac_main <- f
  rm(f, m, t); invisible(gc())
}

market_long <- rbindlist(mkt_list)
themes_long <- rbindlist(th_list)
setorder(market_long, weighting, country, date)
setorder(themes_long, weighting, country, theme, date)

# Sign convention check: directions in the data match the official list?
dir_check <- merge(unique(fac_main[, .(factor, direction)]),
                   factor_dictionary[, .(factor, direction_xlsx)], by = "factor")
n_diff <- dir_check[as.numeric(direction) != direction_xlsx, .N]
say("Factors whose data 'direction' differs from factor_details.xlsx: ", n_diff,
    if (n_diff > 0) paste0(" (", paste(dir_check[as.numeric(direction) != direction_xlsx, factor],
                                        collapse = ", "), ")") else "")

# -----------------------------------------------------------------------------
# 3. Verify how theme returns are built from factor returns
#    JKP documentation: theme = equal-weighted average of the factors in the
#    cluster, after signing each factor by its original paper's direction.
#    We test whether the factor 'ret' column is already signed.
# -----------------------------------------------------------------------------
recon <- fac_main[!is.na(ret), .(
  avg_as_is   = mean(ret),
  avg_flipped = mean(as.numeric(direction) * ret),
  n_used      = .N
), by = .(country, date, theme)]
recon <- merge(recon, th_list[[w_main]][, .(country, date, theme, theme_ret = ret, n_factors)],
               by = c("country", "date", "theme"))

recon_summary <- recon[, .(
  months            = .N,
  corr_as_is        = cor(avg_as_is, theme_ret),
  max_abs_diff_as_is   = max(abs(avg_as_is - theme_ret)),
  corr_flipped      = cor(avg_flipped, theme_ret),
  max_abs_diff_flipped = max(abs(avg_flipped - theme_ret)),
  share_same_n      = mean(n_used == n_factors)
), by = theme]
fwrite(recon_summary, "output/check_theme_reconstruction.csv")

as_is_ok <- recon_summary[, all(max_abs_diff_as_is < 1e-8)]
say("\nTheme reconstruction check (", w_main, "):")
say("  theme = plain average of factor 'ret' column? ", as_is_ok)
say("  largest gap, plain average:   ", signif(recon_summary[, max(max_abs_diff_as_is)], 3))
say("  largest gap, direction*ret:   ", signif(recon_summary[, max(max_abs_diff_flipped)], 3))
say("  -> factor returns in the files are ",
    if (as_is_ok) "ALREADY SIGNED (long = higher expected return side); use 'ret' as is."
    else "NOT a plain average; inspect output/check_theme_reconstruction.csv")

# -----------------------------------------------------------------------------
# 4. Wide panel: one row per country-month, market + 13 themes
# -----------------------------------------------------------------------------
build_panel <- function(w) {
  tw <- dcast(th_list[[w]], country + date ~ theme, value.var = "ret")
  p  <- merge(mkt_list[[w]][, .(country, date, mkt, mkt_n_stocks)], tw,
              by = c("country", "date"), all = TRUE)
  p[, n_themes := rowSums(!is.na(.SD)), .SDcols = themes]
  p[, complete := !is.na(mkt) & n_themes == length(themes)]
  p <- merge(p, countries[, .(country = iso3, country_name = country, msci)],
             by = "country", all.x = TRUE)
  setcolorder(p, c("country", "country_name", "msci", "date", "mkt", "mkt_n_stocks", themes))
  setorder(p, country, date)
  p[]
}

panels <- lapply(setNames(weightings, weightings), build_panel)
for (w in weightings) {
  fwrite(panels[[w]], sprintf("data/panel_%s.csv", w))
  saveRDS(panels[[w]], sprintf("data/panel_%s.rds", w))
}
fwrite(market_long, "data/market_long.csv")
fwrite(themes_long, "data/themes_long.csv")
fwrite(factor_dictionary, "data/factor_dictionary.csv")

P <- panels[[w_main]]
say("\nPanel (", w_main, "): ", nrow(P), " country-months, ", uniqueN(P$country), " countries")
say("  rows with themes but no market return: ", P[is.na(mkt), .N])
say("  rows with market but no theme at all:   ", P[!is.na(mkt) & n_themes == 0, .N])
say("  countries missing from country_meta:    ",
    paste(setdiff(unique(P$country), countries$iso3), collapse = ", "))

# Month gaps inside each country's market series
gap_check <- P[!is.na(mkt), .(gaps = sum(diff(year(date) * 12 + month(date)) != 1)), by = country]
say("  countries with gaps in the market series: ", gap_check[gaps > 0, .N],
    if (gap_check[gaps > 0, .N] > 0) paste0(" (", paste(gap_check[gaps > 0, country], collapse = ", "), ")") else "")

# Units sanity check (returns are decimals: 0.01 = 1%)
us <- P[country == "USA" & !is.na(mkt)]
say(sprintf("  USA market: %s to %s, mean %.2f%%/month, sd %.2f%%/month, annualized Sharpe %.2f",
            min(us$date), max(us$date), 100 * mean(us$mkt), 100 * sd(us$mkt),
            sqrt(12) * mean(us$mkt) / sd(us$mkt)))

# -----------------------------------------------------------------------------
# 5. Coverage tables (inputs to Table 1)
# -----------------------------------------------------------------------------
longest_run <- function(flag) {        # longest stretch of consecutive TRUEs
  r <- rle(flag); if (!any(r$values)) 0L else max(r$lengths[r$values])
}

coverage <- P[, {
  has_mkt <- !is.na(mkt)
  .(country_name     = country_name[1],
    msci             = msci[1],
    mkt_start        = if (any(has_mkt)) min(date[has_mkt]) else as.IDate(NA),
    mkt_end          = if (any(has_mkt)) max(date[has_mkt]) else as.IDate(NA),
    mkt_months       = sum(has_mkt),
    themes_ever      = sum(sapply(.SD, function(x) any(!is.na(x)))),
    complete_start   = if (any(complete)) min(date[complete]) else as.IDate(NA),
    complete_end     = if (any(complete)) max(date[complete]) else as.IDate(NA),
    complete_months  = sum(complete),
    longest_complete_run = longest_run(complete),
    median_stocks    = as.numeric(median(mkt_n_stocks, na.rm = TRUE)),
    median_stocks_last5y = as.numeric(median(mkt_n_stocks[date > max(date) - 5 * 365], na.rm = TRUE)))
}, by = country, .SDcols = themes]
coverage[, eligible := complete_months >= min_months]
setorder(coverage, -complete_months)
fwrite(coverage, "output/coverage_by_country.csv")

theme_cov <- dcast(th_list[[w_main]][!is.na(ret), .N, by = .(country, theme)],
                   country ~ theme, value.var = "N", fill = 0)
fwrite(theme_cov, "output/coverage_theme_by_country.csv")

say("\nCoverage: ", coverage[eligible == TRUE, .N], " of ", nrow(coverage),
    " countries have >= ", min_months, " complete months (market + all 13 themes)")
say("  by MSCI group (eligible): ",
    paste(coverage[eligible == TRUE, .N, by = msci][, paste0(msci, " ", N)], collapse = ", "))

# -----------------------------------------------------------------------------
# 6. Common-history window for the five pilot countries
# -----------------------------------------------------------------------------
pp <- P[country %in% pilot & complete == TRUE]
common_dates <- Reduce(intersect, split(pp$date, pp$country))
common_dates <- sort(as.IDate(common_dates))
cw <- data.table(
  countries     = paste(pilot, collapse = " "),
  common_start  = min(common_dates),
  common_end    = max(common_dates),
  common_months = length(common_dates),
  contiguous    = length(common_dates) ==
                  (year(max(common_dates)) * 12 + month(max(common_dates))) -
                  (year(min(common_dates)) * 12 + month(min(common_dates))) + 1
)
fwrite(cw, "output/common_window_feasibility.csv")
say("\nCommon window, pilot countries: ", cw$common_start, " to ", cw$common_end,
    " (", cw$common_months, " months, contiguous = ", cw$contiguous, ")")
print(coverage[country %in% pilot,
               .(country, mkt_start, complete_start, complete_end, complete_months, median_stocks)])

say("\nDone. Files written to data/ and output/.")
