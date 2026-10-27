# =============================================================================
# 05_importance.R  --  Model 2: relative explanatory importance (Tables 4 and 5)
# FIN 864 project: International Heterogeneity in Equity Factor Exposures
#
# Run from the project root (after 02_build_panel.R):
#   Rscript scripts/05_importance.R
#
# For each pilot country, over the common-history window:
#   (a) standardized coefficients (Y and X scaled to unit variance)
#   (b) drop-one incremental R^2:  dR2_k = R2_full - R2_without_k          (Table 4)
#   (c) exact Shapley decomposition of R^2 over all 2^13 theme subsets   (Table 5)
#       reported as percentage points of total variance and as % of explained
#   (d) stability: moving-block bootstrap of Shapley shares and rankings
#
# R^2 of any subset S is computed from correlations:  R2(S) = r_S' R_SS^-1 r_S,
# which equals the OLS R^2 with an intercept.
# =============================================================================

suppressPackageStartupMessages(library(data.table))
set.seed(864)

pilot  <- c("USA", "JPN", "DEU", "GBR", "AUS")
themes <- c("accruals", "debt_issuance", "investment", "low_leverage", "low_risk",
            "momentum", "profit_growth", "profitability", "quality", "seasonality",
            "short_term_reversal", "size", "value")
theme_label <- c(accruals = "Accruals", debt_issuance = "Debt issuance", investment = "Investment",
                 low_leverage = "Low leverage", low_risk = "Low risk", momentum = "Momentum",
                 profit_growth = "Profit growth", profitability = "Profitability", quality = "Quality",
                 seasonality = "Seasonality", short_term_reversal = "Short-term reversal",
                 size = "Size", value = "Value")
K <- length(themes)
B_REPS    <- 200   # bootstrap replications
BLOCK_LEN <- 12    # months per block

P  <- readRDS("data/panel_vw_cap.rds")
pp <- P[country %in% pilot & complete == TRUE]
common_dates <- sort(as.IDate(Reduce(intersect, split(pp$date, pp$country))))

# ---- all subsets of the 13 themes, encoded as bitmasks 0 .. 2^K - 1 ----------
n_sub   <- 2^K
members <- lapply(0:(n_sub - 1), function(m) which(bitwAnd(m, 2^(0:(K - 1))) > 0))
sizes   <- lengths(members)

subset_r2 <- function(Rxx, rxy) {
  out <- numeric(n_sub)
  for (m in 2:n_sub) {                      # m = 1 is the empty set, R2 = 0
    s <- members[[m]]
    out[m] <- sum(rxy[s] * solve(Rxx[s, s, drop = FALSE], rxy[s]))
  }
  out
}

# Shapley weights: |S|! (K - |S| - 1)! / K!  for S not containing k
w_size <- factorial(0:(K - 1)) * factorial((K - 1):0) / factorial(K)
shapley_from_r2 <- function(r2) {
  sapply(seq_len(K), function(k) {
    bit <- 2^(k - 1)
    without <- which(bitwAnd(0:(n_sub - 1), bit) == 0)      # indices (1-based) of S without k
    with    <- without + bit                                 # same S plus k
    sum(w_size[sizes[without] + 1] * (r2[with] - r2[without]))
  })
}

importance <- function(d) {
  X <- as.matrix(d[, ..themes]); y <- d$mkt
  R <- cor(cbind(y, X)); Rxx <- R[-1, -1]; rxy <- R[-1, 1]
  r2  <- subset_r2(Rxx, rxy)
  full <- r2[n_sub]
  drop_one <- sapply(seq_len(K), function(k) full - r2[n_sub - 2^(k - 1)])
  shap <- shapley_from_r2(r2)
  std_beta <- solve(Rxx, rxy)               # standardized OLS coefficients
  list(full = full, drop_one = drop_one, shap = shap, std_beta = std_beta)
}

block_boot_idx <- function(n, L) {          # moving-block bootstrap indices
  starts <- sample.int(n - L + 1, ceiling(n / L), replace = TRUE)
  idx <- as.vector(sapply(starts, function(s) s:(s + L - 1)))
  idx[seq_len(n)]
}

rows <- list(); boot_rows <- list()
for (cty in pilot) {
  d   <- pp[country == cty & date %in% common_dates][order(date)]
  imp <- importance(d)
  stopifnot(abs(sum(imp$shap) - imp$full) < 1e-10)          # Shapley adds up to R2
  rows[[cty]] <- data.table(country = cty, theme = themes, full_r2 = imp$full,
                            std_beta = imp$std_beta, drop_one = imp$drop_one,
                            shapley = imp$shap, shapley_share = imp$shap / imp$full,
                            shapley_rank = rank(-imp$shap))
  message("Bootstrapping ", cty, " ...")
  bs <- t(replicate(B_REPS, {
    i <- block_boot_idx(nrow(d), BLOCK_LEN)
    s <- importance(d[i])$shap
    c(s / sum(s), rank(-s))
  }))
  boot_rows[[cty]] <- data.table(country = cty, theme = themes,
    share_p05 = apply(bs[, 1:K], 2, quantile, 0.05),
    share_p95 = apply(bs[, 1:K], 2, quantile, 0.95),
    rank_p05  = apply(bs[, K + 1:K], 2, quantile, 0.05),
    rank_p95  = apply(bs[, K + 1:K], 2, quantile, 0.95),
    prob_top1 = colMeans(bs[, K + 1:K] == 1),
    prob_top3 = colMeans(bs[, K + 1:K] <= 3))
}
imp_all  <- rbindlist(rows)
boot_all <- rbindlist(boot_rows)
res <- merge(imp_all, boot_all, by = c("country", "theme"))
fwrite(res[, lapply(.SD, function(v) if (is.numeric(v)) signif(v, 5) else v)], "output/table4_5_importance_long.csv")

# ---- formatted tables ---------------------------------------------------------
wide <- function(expr) {
  x <- dcast(res[, .(theme, country, v = eval(expr))], theme ~ country, value.var = "v")
  setcolorder(x, c("theme", pilot)); x[match(themes, theme)][, theme := theme_label[theme]]
}
full_r2 <- res[, .(r2 = full_r2[1]), by = country][match(pilot, country)]
sum_drop <- res[, .(s = sum(drop_one)), by = country][match(pilot, country)]

t4_std  <- wide(quote(sprintf("%.2f", std_beta)))
t4_drop <- wide(quote(sprintf("%.1f", 100 * drop_one)))
t4_drop <- rbind(t4_drop,
  data.table(theme = "Sum of drop-one", t(setNames(sprintf("%.1f", 100 * sum_drop$s), pilot))),
  data.table(theme = "Full-model R2",   t(setNames(sprintf("%.1f", 100 * full_r2$r2), pilot))))
t5_abs   <- wide(quote(sprintf("%.1f", 100 * shapley)))
t5_abs   <- rbind(t5_abs, data.table(theme = "Total = R2", t(setNames(sprintf("%.1f", 100 * full_r2$r2), pilot))))
t5_share <- wide(quote(sprintf("%.0f%% [%.0f-%.0f]", 100 * shapley_share, 100 * share_p05, 100 * share_p95)))
t5_rank  <- wide(quote(sprintf("%d (%.0f%%)", as.integer(shapley_rank), 100 * prob_top3)))

for (nm in c("t4_std", "t4_drop", "t5_abs", "t5_share", "t5_rank")) fwrite(get(nm), sprintf("output/%s.csv", nm))

sink("output/tables_4_5.txt")
cat("MODEL 2. Relative explanatory importance of the 13 themes, by country\n")
cat("Common-history window:", format(min(common_dates), "%b %Y"), "to", format(max(common_dates), "%b %Y"),
    "(", length(common_dates), "months ). Capped value-weighted USD returns.\n\n")
cat("TABLE 4, Panel A. Standardized coefficients (market and themes scaled to unit variance)\n")
print(t4_std, row.names = FALSE)
cat("\nTABLE 4, Panel B. Drop-one incremental R2 (percentage points of total market variance)\n")
cat("R2 lost when the theme is removed from the full model. Overlap between themes is not counted,\n")
cat("so the column sums are below the full-model R2.\n")
print(t4_drop, row.names = FALSE)
cat("\nTABLE 5, Panel A. Shapley R2 (percentage points of total market variance; columns sum to R2)\n")
print(t5_abs, row.names = FALSE)
cat("\nTABLE 5, Panel B. Shapley share of explained variance, with 90% block-bootstrap interval\n")
cat("(", B_REPS, "replications,", BLOCK_LEN, "-month blocks )\n")
print(t5_share, row.names = FALSE)
cat("\nTABLE 5, Panel C. Shapley rank (1 = most important), and bootstrap probability of ranking in the top 3\n")
print(t5_rank, row.names = FALSE)
sink()
cat(readLines("output/tables_4_5.txt"), sep = "\n")
