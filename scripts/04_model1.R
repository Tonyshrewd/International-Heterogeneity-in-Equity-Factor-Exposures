# =============================================================================
# 04_model1.R  --  Model 1: contemporaneous country-level theme exposures (Table 3)
# FIN 864 project: International Heterogeneity in Equity Factor Exposures
#
# Run from the project root (after 02_build_panel.R):
#   Rscript scripts/04_model1.R
#
# For each pilot country c, estimate by OLS
#   MKT_{c,t} = alpha_c + sum_k beta_{k,c} * THEME_{k,c,t} + e_{c,t}
# on (1) the common-history window shared by all five countries and
#    (2) each country's maximal history (all months with market + 13 themes).
#
# Inference: Newey-West (HAC) standard errors, lag = NW_LAG months.
# Diagnostics: adjusted R^2, sample size, variance inflation factors (VIF),
#              condition number, theme correlation matrices.
# Interpretation: associations only. Market and theme portfolios share stocks.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sandwich)
})

pilot  <- c("USA", "JPN", "DEU", "GBR", "AUS")
themes <- c("accruals", "debt_issuance", "investment", "low_leverage", "low_risk",
            "momentum", "profit_growth", "profitability", "quality", "seasonality",
            "short_term_reversal", "size", "value")
theme_label <- c(accruals = "Accruals", debt_issuance = "Debt issuance", investment = "Investment",
                 low_leverage = "Low leverage", low_risk = "Low risk", momentum = "Momentum",
                 profit_growth = "Profit growth", profitability = "Profitability", quality = "Quality",
                 seasonality = "Seasonality", short_term_reversal = "Short-term reversal",
                 size = "Size", value = "Value")
NW_LAG <- 6   # about 4*(T/100)^(2/9) for T = 419; robustness to other lags below

P  <- readRDS("data/panel_vw_cap.rds")
pp <- P[country %in% pilot & complete == TRUE]
common_dates <- sort(as.IDate(Reduce(intersect, split(pp$date, pp$country))))

fmla <- as.formula(paste("mkt ~", paste(themes, collapse = " + ")))

vif <- function(X) {                       # VIF_j = 1 / (1 - R^2 of x_j on other x's)
  sapply(colnames(X), function(j) {
    r2 <- summary(lm(X[, j] ~ X[, setdiff(colnames(X), j)]))$r.squared
    1 / (1 - r2)
  })
}

fit_country <- function(d, cty, sample_name) {
  fit <- lm(fmla, data = d)
  V   <- NeweyWest(fit, lag = NW_LAG, prewhite = FALSE, adjust = TRUE)
  se  <- sqrt(diag(V))
  b   <- coef(fit)
  X   <- as.matrix(d[, ..themes])
  v   <- vif(X)
  ev  <- eigen(cor(X), only.values = TRUE)$values
  coefs <- data.table(country = cty, sample = sample_name, term = names(b),
                      beta = b, se_nw = se, t_nw = b / se,
                      p_nw = 2 * pt(-abs(b / se), df = fit$df.residual),
                      vif = c(NA, v[names(b)[-1]]))
  # Robustness of inference to the NW lag choice
  lag_t <- sapply(c(3, 12), function(L) b / sqrt(diag(NeweyWest(fit, lag = L, prewhite = FALSE, adjust = TRUE))))
  coefs[, `:=`(t_nw_lag3 = lag_t[, 1], t_nw_lag12 = lag_t[, 2])]
  diag_row <- data.table(country = cty, sample = sample_name,
                         start = min(d$date), end = max(d$date), N = nobs(fit),
                         r2 = summary(fit)$r.squared, adj_r2 = summary(fit)$adj.r.squared,
                         max_vif = max(v), max_vif_theme = names(which.max(v)),
                         mean_vif = mean(v), condition_number = sqrt(max(ev) / min(ev)))
  list(coefs = coefs, diag = diag_row, cor = cor(X))
}

res <- list()
for (cty in pilot) {
  d_full <- pp[country == cty][order(date)]
  res[[paste(cty, "common")]]  <- fit_country(d_full[date %in% common_dates], cty, "common")
  res[[paste(cty, "maximal")]] <- fit_country(d_full, cty, "maximal")
}
coefs <- rbindlist(lapply(res, `[[`, "coefs"))
diags <- rbindlist(lapply(res, `[[`, "diag"))
fwrite(coefs[, lapply(.SD, function(v) if (is.numeric(v)) signif(v, 5) else v)], "output/table3_model1_coefficients_long.csv")
fwrite(diags[, lapply(.SD, function(v) if (is.numeric(v)) signif(v, 5) else v)], "output/table3_model1_diagnostics.csv")

# Theme correlation matrices (diagnostic memo)
for (cty in pilot) {
  cm <- round(res[[paste(cty, "common")]]$cor, 2)
  fwrite(data.table(theme = rownames(cm), cm), sprintf("output/theme_correlations_%s.csv", cty))
}

# -----------------------------------------------------------------------------
# Formatted Table 3 (theme rows, country columns): beta (NW t-stat) and stars
# -----------------------------------------------------------------------------
stars <- function(p) fifelse(p < 0.01, "***", fifelse(p < 0.05, "**", fifelse(p < 0.10, "*", "")))
make_table <- function(s) {
  x <- coefs[sample == s & term != "(Intercept)",
             .(term, country, cell = sprintf("%.2f%s (%.1f)", beta, stars(p_nw), t_nw))]
  w <- dcast(x, term ~ country, value.var = "cell")
  setcolorder(w, c("term", pilot))
  w <- w[match(themes, term)][, term := theme_label[term]]
  dg <- diags[sample == s][match(pilot, country)]
  extra <- rbind(
    data.table(term = "Adjusted R2", t(setNames(sprintf("%.2f", dg$adj_r2), pilot))),
    data.table(term = "Months", t(setNames(as.character(dg$N), pilot))),
    data.table(term = "Max VIF", t(setNames(sprintf("%.1f", dg$max_vif), pilot))),
    data.table(term = "Condition number", t(setNames(sprintf("%.1f", dg$condition_number), pilot))))
  rbind(w, extra, use.names = TRUE)
}
t3_common  <- make_table("common")
t3_maximal <- make_table("maximal")
fwrite(t3_common,  "output/table3_model1_common.csv")
fwrite(t3_maximal, "output/table3_model1_maximal.csv")

sink("output/table3_model1.txt")
cat("TABLE 3. Model 1: market excess return on 13 theme returns, country by country\n")
cat("Cells: OLS coefficient, stars from Newey-West p-values (lag", NW_LAG, "), NW t-statistic in parentheses\n")
cat("*** p<0.01, ** p<0.05, * p<0.10. Capped value-weighted USD returns.\n\n")
cat("Panel A. Common-history window,", format(min(common_dates), "%b %Y"), "to",
    format(max(common_dates), "%b %Y"), "\n")
print(t3_common, row.names = FALSE)
cat("\nPanel B. Maximal history for each country\n")
print(t3_maximal, row.names = FALSE)
cat("\nSample windows (maximal):\n")
print(diags[sample == "maximal", .(country, start, end, N)], row.names = FALSE)
cat("\nLargest VIF per country (common window):\n")
print(diags[sample == "common", .(country, max_vif = round(max_vif, 1), max_vif_theme, mean_vif = round(mean_vif, 1))], row.names = FALSE)
cat("\nSign changes in t-statistics across NW lags 3 / 6 / 12 (common window): ",
    coefs[sample == "common" & term != "(Intercept)",
          sum((abs(t_nw) > 1.96) != (abs(t_nw_lag3) > 1.96) | (abs(t_nw) > 1.96) != (abs(t_nw_lag12) > 1.96))],
    "coefficients change 5% significance status\n")
sink()
cat(readLines("output/table3_model1.txt"), sep = "\n")
