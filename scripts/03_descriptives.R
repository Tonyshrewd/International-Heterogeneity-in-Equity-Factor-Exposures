# =============================================================================
# 03_descriptives.R  --  Tables 1 and 2 for the five-country feasibility study
# FIN 864 project: International Heterogeneity in Equity Factor Exposures
#
# Run from the project root (after 02_build_panel.R):
#   Rscript scripts/03_descriptives.R
#
# Table 1  Sample construction: coverage, windows, number of stocks
# Table 2  Descriptive statistics over the common-history window
#          Panel A  market excess returns
#          Panel B  average theme returns (premia), t-statistics, Sharpe ratios
#          Panel C  correlation of each theme with the country market return
#
# All returns are monthly, in USD, capped value weighted (JKP "vw_cap").
# Outputs: output/table1_*.csv, output/table2*_*.csv, output/tables_1_2.txt
# =============================================================================

suppressPackageStartupMessages(library(data.table))

pilot  <- c("USA", "JPN", "DEU", "GBR", "AUS")
themes <- c("accruals", "debt_issuance", "investment", "low_leverage", "low_risk",
            "momentum", "profit_growth", "profitability", "quality", "seasonality",
            "short_term_reversal", "size", "value")
theme_label <- c(accruals = "Accruals", debt_issuance = "Debt issuance", investment = "Investment",
                 low_leverage = "Low leverage", low_risk = "Low risk", momentum = "Momentum",
                 profit_growth = "Profit growth", profitability = "Profitability", quality = "Quality",
                 seasonality = "Seasonality", short_term_reversal = "Short-term reversal",
                 size = "Size", value = "Value")

P <- readRDS("data/panel_vw_cap.rds")
pp <- P[country %in% pilot]

# Common-history window: months where all five countries have the market and all 13 themes
common_dates <- sort(as.IDate(Reduce(intersect, split(pp[complete == TRUE]$date,
                                                      pp[complete == TRUE]$country))))
W <- pp[date %in% common_dates]
stopifnot(W[, all(complete)], W[, .N, by = country][, all(N == length(common_dates))])

ym <- function(d) format(d, "%b %Y")

# -----------------------------------------------------------------------------
# Table 1: sample construction
# -----------------------------------------------------------------------------
table1 <- pp[, {
  has_mkt <- !is.na(mkt)
  .(Country            = country_name[1],
    MSCI               = msci[1],
    `Market from`      = ym(min(date[has_mkt])),
    `All 13 themes from` = ym(min(date[complete])),
    `Last month`       = ym(max(date[complete])),
    `Months (maximal)` = sum(complete),
    `Months (common)`  = sum(date %in% common_dates),
    `Median stocks`    = as.integer(median(mkt_n_stocks[date %in% common_dates])),
    `Min stocks`       = as.integer(min(mkt_n_stocks[date %in% common_dates])))
}, by = country][match(pilot, country)]
setnames(table1, "country", "Code")
fwrite(table1, "output/table1_sample.csv")

# Appendix: all countries with at least 120 complete months
cov <- fread("output/coverage_by_country.csv")
fwrite(cov[eligible == TRUE, .(country, country_name, msci, complete_start, complete_end,
                               complete_months, longest_complete_run, median_stocks)],
       "output/table1_appendix_eligible_countries.csv")

# -----------------------------------------------------------------------------
# Table 2: descriptive statistics over the common window
# -----------------------------------------------------------------------------
stats <- function(x) {
  n <- length(x); m <- mean(x); s <- sd(x)
  list(n = n, mean = 100 * m, sd = 100 * s, t = m / (s / sqrt(n)),
       sharpe = sqrt(12) * m / s, min = 100 * min(x), max = 100 * max(x),
       skew = mean((x - m)^3) / s^3, ac1 = cor(x[-1], x[-n]))
}

# Panel A: market excess return
t2a <- W[, {
  s <- stats(mkt)
  .(Country = country_name[1], Months = s$n,
    `Mean (%/mo)` = round(s$mean, 2), `SD (%/mo)` = round(s$sd, 2), `t(mean)` = round(s$t, 2),
    `Ann. Sharpe` = round(s$sharpe, 2), `Min (%)` = round(s$min, 1), `Max (%)` = round(s$max, 1),
    Skewness = round(s$skew, 2), `AR(1)` = round(s$ac1, 2))
}, by = country][match(pilot, country)]
fwrite(t2a, "output/table2a_market.csv")

# Panel B: average theme returns (premia), long format and wide tables
long <- melt(W, id.vars = c("country", "date", "mkt"), measure.vars = themes,
             variable.name = "theme", value.name = "ret")
long[, theme := as.character(theme)]
prem <- long[, {
  s <- stats(ret)
  .(mean = s$mean, sd = s$sd, t = s$t, sharpe = s$sharpe, corr_mkt = cor(ret, mkt))
}, by = .(country, theme)]

wide_fmt <- function(dt, value_expr) {
  x <- dcast(dt[, .(theme, country, v = eval(value_expr))], theme ~ country, value.var = "v")
  setcolorder(x, c("theme", pilot))
  x[, theme := theme_label[theme]][order(match(theme, theme_label))]
}
t2b_mean   <- wide_fmt(prem, quote(sprintf("%.2f (%.2f)", mean, t)))
t2b_sharpe <- wide_fmt(prem, quote(round(sharpe, 2)))
t2c_corr   <- wide_fmt(prem, quote(round(corr_mkt, 2)))
fwrite(prem[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 4) else v)],
       "output/table2_theme_stats_long.csv")
fwrite(t2b_mean,   "output/table2b_theme_premia.csv")
fwrite(t2b_sharpe, "output/table2b_theme_sharpe.csv")
fwrite(t2c_corr,   "output/table2c_theme_market_corr.csv")

# -----------------------------------------------------------------------------
# Plain-text version of both tables
# -----------------------------------------------------------------------------
sink("output/tables_1_2.txt")
cat("Five-country feasibility study, JKP Global Factor Data (capped value weighted, USD)\n")
cat("Common-history window:", ym(min(common_dates)), "to", ym(max(common_dates)),
    "(", length(common_dates), "months )\n\n")
cat("TABLE 1. Sample construction\n"); print(table1, row.names = FALSE)
cat("\nTABLE 2, Panel A. Market excess returns (common window)\n"); print(t2a, row.names = FALSE)
cat("\nTABLE 2, Panel B. Average theme returns, % per month (t-statistic), common window\n")
print(t2b_mean, row.names = FALSE)
cat("\nTABLE 2, Panel B (cont.). Annualized Sharpe ratios of theme returns\n")
print(t2b_sharpe, row.names = FALSE)
cat("\nTABLE 2, Panel C. Correlation of each theme return with the country market return\n")
print(t2c_corr, row.names = FALSE)
sink()

cat(readLines("output/tables_1_2.txt"), sep = "\n")
