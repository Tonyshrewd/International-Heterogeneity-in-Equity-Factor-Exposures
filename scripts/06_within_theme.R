# =============================================================================
# 06_within_theme.R  --  Model 3: individual factors within themes (Table 6)
# FIN 864 project: International Heterogeneity in Equity Factor Exposures
#
# Run from the project root (after 02_build_panel.R and 05_importance.R):
#   Rscript scripts/06_within_theme.R
#
# For the most important themes (Low risk in every country, plus each
# country's second-ranked theme by Shapley R2 from Table 5), regress the
# country market excess return on the theme's individual factors:
#   MKT_{c,t} = a + sum_{j in J_k} gamma_{j,c} F_{j,c,t} + u_{c,t}
# and decompose the R2 across factors with an exact Shapley decomposition.
#
# Sample rules (common-history window, Feb 1991 - Dec 2025):
#   * keep a factor if it is observed in >= 85% of window months;
#   * for Low risk, use the factors that pass in ALL five countries, so the
#     factor rankings are comparable across countries (H4);
#   * estimate on months where every kept factor is observed.
# Also reported: R2 of the theme average alone, and the share of factor
# variance captured by the first principal component (how one-dimensional
# the theme is). Results are descriptive: factors within a theme overlap.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

pilot <- c("USA", "JPN", "DEU", "GBR", "AUS")
MIN_COVERAGE <- 0.85
theme_label <- c(accruals = "Accruals", debt_issuance = "Debt issuance", investment = "Investment",
                 low_leverage = "Low leverage", low_risk = "Low risk", momentum = "Momentum",
                 profit_growth = "Profit growth", profitability = "Profitability", quality = "Quality",
                 seasonality = "Seasonality", short_term_reversal = "Short-term reversal",
                 size = "Size", value = "Value")

P    <- readRDS("data/panel_vw_cap.rds")
F    <- readRDS("data/factors_long_vw_cap.rds")
dict <- fread("data/factor_dictionary.csv")
pp   <- P[country %in% pilot & complete == TRUE]
common_dates <- sort(as.IDate(Reduce(intersect, split(pp$date, pp$country))))
n_window <- length(common_dates)

# Theme ranking from Table 5 decides which themes to open up
imp <- fread("output/table4_5_importance_long.csv")
second <- imp[shapley_rank == 2, .(country, theme)]
plan <- rbind(data.table(country = pilot, theme = "low_risk"), second)[order(match(country, pilot))]

# ---- exact Shapley decomposition for K factors (generic K) --------------------
shapley_exact <- function(Rxx, rxy) {
  K <- length(rxy); n_sub <- 2^K
  bits <- 2^(0:(K - 1))
  r2 <- numeric(n_sub); size <- integer(n_sub)
  for (m in 1:(n_sub - 1)) {                 # mask m (0-based), stored at m + 1
    s <- which(bitwAnd(m, bits) > 0)
    size[m + 1] <- length(s)
    r2[m + 1] <- sum(rxy[s] * solve(Rxx[s, s, drop = FALSE], rxy[s]))
  }
  w <- factorial(0:(K - 1)) * factorial((K - 1):0) / factorial(K)
  masks <- 0:(n_sub - 1)
  phi <- sapply(seq_len(K), function(k) {
    without <- which(bitwAnd(masks, bits[k]) == 0)
    sum(w[size[without] + 1] * (r2[without + bits[k]] - r2[without]))
  })
  list(phi = phi, r2_full = r2[n_sub])
}

kept_factors <- function(cty, th) {
  x <- F[country == cty & theme == th & date %in% common_dates & !is.na(ret), .N, by = factor]
  x[N >= MIN_COVERAGE * n_window, factor]
}
# Low risk: common factor set across the five countries
lr_common <- Reduce(intersect, lapply(pilot, kept_factors, th = "low_risk"))

summary_rows <- list(); factor_rows <- list()
for (i in seq_len(nrow(plan))) {
  cty <- plan$country[i]; th <- plan$theme[i]
  all_f  <- dict[theme == th, factor]
  keep   <- if (th == "low_risk") lr_common else kept_factors(cty, th)
  wide   <- dcast(F[country == cty & theme == th & date %in% common_dates & factor %in% keep],
                  date ~ factor, value.var = "ret")
  d <- merge(pp[country == cty, .(date, mkt, theme_ret = get(th))], wide, by = "date")
  d <- d[complete.cases(d)]
  X <- as.matrix(d[, ..keep]); y <- d$mkt
  R <- cor(cbind(y, X)); Rxx <- R[-1, -1]; rxy <- R[-1, 1]
  message(sprintf("%s %s: %d factors, %d months", cty, th, length(keep), nrow(d)))
  sh <- shapley_exact(Rxx, rxy)
  stopifnot(abs(sum(sh$phi) - sh$r2_full) < 1e-8)
  ev <- eigen(Rxx, only.values = TRUE)$values

  summary_rows[[i]] <- data.table(
    country = cty, theme = th, rank_in_table5 = if (th == "low_risk") 1L else 2L,
    factors_available = length(all_f), factors_used = length(keep), months = nrow(d),
    first_month = min(d$date), r2_theme_average = cor(y, d$theme_ret)^2,
    r2_all_factors = sh$r2_full, pc1_share = ev[1] / sum(ev),
    dropped = paste(setdiff(all_f, keep), collapse = " "))
  factor_rows[[i]] <- data.table(country = cty, theme = th, factor = keep,
                                 corr_mkt = rxy, shapley = sh$phi,
                                 share = sh$phi / sh$r2_full)[order(-shapley)][, rank := .I]
}
summ <- rbindlist(summary_rows)
facs <- rbindlist(factor_rows)
facs <- merge(facs, dict[, .(factor, description)], by = "factor", all.x = TRUE)[order(match(country, pilot), theme, rank)]
fwrite(summ[, lapply(.SD, function(v) if (is.numeric(v)) signif(v, 4) else v)], "output/table6_within_theme_summary.csv")
fwrite(facs[, lapply(.SD, function(v) if (is.numeric(v)) signif(v, 4) else v)], "output/table6_within_theme_factors.csv")

# ---- formatted output -----------------------------------------------------------
top3 <- facs[rank <= 3, .(top = paste(sprintf("%s (%.0f%%)", description, 100 * share), collapse = "; ")),
             by = .(country, theme)]
t6 <- merge(summ, top3, by = c("country", "theme"))[order(match(country, pilot), -rank_in_table5)]
t6_print <- t6[, .(Country = country, Theme = theme_label[theme], `Theme rank` = rank_in_table5,
                   `Factors used/avail.` = paste0(factors_used, "/", factors_available), Months = months,
                   `R2 theme avg` = sprintf("%.2f", r2_theme_average),
                   `R2 all factors` = sprintf("%.2f", r2_all_factors),
                   `PC1 share` = sprintf("%.2f", pc1_share), `Top 3 factors (Shapley share)` = top)]
fwrite(t6_print, "output/table6_within_theme.csv")

lr_rank <- dcast(facs[theme == "low_risk"], description ~ country, value.var = "rank")
setcolorder(lr_rank, c("description", pilot))
lr_rank[, avg := rowMeans(.SD), .SDcols = pilot]
lr_rank <- lr_rank[order(avg)][, avg := NULL]
fwrite(lr_rank, "output/table6_low_risk_factor_ranks.csv")

sink("output/table6_within_theme.txt")
cat("TABLE 6. Model 3: individual factors within the most important themes\n")
cat("Market excess return regressed on the theme's individual factors; exact Shapley decomposition of R2.\n")
cat("Common window", format(min(common_dates), "%b %Y"), "to", format(max(common_dates), "%b %Y"),
    "; factors kept if observed in >=", 100 * MIN_COVERAGE, "% of months; complete months only.\n")
cat("Low risk uses the same", length(lr_common), "factors in every country.\n\n")
for (r in seq_len(nrow(t6_print))) {
  x <- t6_print[r]
  cat(sprintf("%s | %s (theme rank %d) | factors %s | %d months\n", x$Country, x$Theme, x$`Theme rank`,
              x$`Factors used/avail.`, x$Months))
  cat(sprintf("   R2: theme average %s, all factors %s | first principal component explains %s of factor variance\n",
              x$`R2 theme avg`, x$`R2 all factors`, x$`PC1 share`))
  cat(sprintf("   Top factors: %s\n", x$`Top 3 factors (Shapley share)`))
  if (nchar(t6[r]$dropped) > 0) cat(sprintf("   Not used (coverage): %s\n", t6[r]$dropped))
  cat("\n")
}
cat("Low risk: rank of each factor by Shapley R2 within the theme (1 = most important)\n")
print(lr_rank, row.names = FALSE)
sink()
cat(readLines("output/table6_within_theme.txt"), sep = "\n")
