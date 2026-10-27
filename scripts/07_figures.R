# =============================================================================
# 07_figures.R  --  Figure 1: theme importance across countries (heatmaps)
# FIN 864 project: International Heterogeneity in Equity Factor Exposures
#
# Run from the project root (after 05_importance.R):
#   Rscript scripts/07_figures.R
#
# Figure 1A  Shapley share of explained market variance, 13 themes x 5 countries
#            (from output/table4_5_importance_long.csv)
# Figure 1B  Same decomposition for the other 12 themes after controlling for
#            Low risk: the market and the 12 themes are residualized on the Low
#            risk theme, and the residual R2 is split by Shapley.
#            This shows the cross-country differences that Low risk otherwise swamps.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

pilot  <- c("USA", "JPN", "DEU", "GBR", "AUS")
cty_label <- c(USA = "United States", JPN = "Japan", DEU = "Germany", GBR = "United Kingdom", AUS = "Australia")
themes <- c("accruals", "debt_issuance", "investment", "low_leverage", "low_risk",
            "momentum", "profit_growth", "profitability", "quality", "seasonality",
            "short_term_reversal", "size", "value")
theme_label <- c(accruals = "Accruals", debt_issuance = "Debt issuance", investment = "Investment",
                 low_leverage = "Low leverage", low_risk = "Low risk", momentum = "Momentum",
                 profit_growth = "Profit growth", profitability = "Profitability", quality = "Quality",
                 seasonality = "Seasonality", short_term_reversal = "Short-term reversal",
                 size = "Size", value = "Value")
others <- setdiff(themes, "low_risk")

# Sequential blue ramp (light = near zero, dark = large)
ramp <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
INK <- "#1F2328"; MUTED <- "#59636E"

# ---- Panel B: Shapley among 12 themes after controlling for Low risk ----------
P  <- readRDS("data/panel_vw_cap.rds")
pp <- P[country %in% pilot & complete == TRUE]
common_dates <- sort(as.IDate(Reduce(intersect, split(pp$date, pp$country))))

shapley_exact <- function(Rxx, rxy) {
  K <- length(rxy); n_sub <- 2^K; bits <- 2^(0:(K - 1))
  r2 <- numeric(n_sub); size <- integer(n_sub)
  for (m in 1:(n_sub - 1)) {
    s <- which(bitwAnd(m, bits) > 0); size[m + 1] <- length(s)
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

partial <- rbindlist(lapply(pilot, function(cty) {
  d <- pp[country == cty & date %in% common_dates]
  resid_on_lr <- function(v) resid(lm(v ~ d$low_risk))
  y <- resid_on_lr(d$mkt)
  X <- sapply(others, function(th) resid_on_lr(d[[th]]))
  R <- cor(cbind(y, X)); sh <- shapley_exact(R[-1, -1], R[-1, 1])
  data.table(country = cty, theme = others, shapley = sh$phi,
             share = sh$phi / sh$r2_full, partial_r2 = sh$r2_full)
}))
fwrite(partial[, lapply(.SD, function(v) if (is.numeric(v)) signif(v, 4) else v)],
       "output/figure1b_shapley_after_low_risk.csv")

# ---- plotting -------------------------------------------------------------------
heat <- function(dt, title, subtitle, legend_title, file, caption) {
  ord <- dt[, .(m = mean(share)), by = theme][order(m), theme]          # most important on top
  dt  <- copy(dt)[, `:=`(theme_f = factor(theme_label[theme], levels = theme_label[ord]),
                         cty_f = factor(cty_label[country], levels = cty_label[pilot]))]
  lim <- c(0, ceiling(max(dt$share) * 20) / 20)
  dt[, ink := fifelse(share > 0.55 * lim[2], "#FFFFFF", INK)]
  g <- ggplot(dt, aes(cty_f, theme_f, fill = share)) +
    geom_tile(color = "#FFFFFF", linewidth = 1.2) +
    geom_text(aes(label = sprintf("%.0f%%", 100 * share), color = ink), size = 3.4) +
    scale_color_identity() +
    scale_fill_gradientn(colours = ramp, limits = lim, labels = function(x) paste0(round(100 * x), "%"),
                         name = legend_title) +
    scale_x_discrete(position = "top") +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL, caption = caption) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", color = INK, size = 13),
          plot.subtitle = element_text(color = MUTED, size = 10, margin = margin(b = 8)),
          plot.caption = element_text(color = MUTED, size = 8, hjust = 0),
          axis.text = element_text(color = INK), axis.text.x.top = element_text(face = "bold"),
          legend.title = element_text(color = MUTED, size = 9), legend.text = element_text(color = MUTED, size = 8),
          plot.background = element_rect(fill = "#FFFFFF", color = NA),
          plot.title.position = "plot", plot.caption.position = "plot",
          plot.margin = margin(12, 14, 10, 10))
  ggsave(file, g, width = 7.2, height = 6.2, dpi = 300, bg = "#FFFFFF")
}

imp <- fread("output/table4_5_importance_long.csv")[, .(country, theme, share = shapley_share)]
heat(imp,
     "Figure 1A. Which themes explain national market returns?",
     "Shapley share of the explained variance of each country's\nmarket excess return across the 13 themes",
     "Share of\nexplained\nvariance", "output/figure1a_theme_importance.png",
     "Monthly capped value-weighted USD returns, Feb 1991-Dec 2025 (419 months).\nEach column sums to 100%.")

r2b <- partial[, .(r2 = partial_r2[1]), by = country][match(pilot, country)]
heat(partial[, .(country, theme, share)],
     "Figure 1B. The other themes, after controlling for Low risk",
     "Shapley share of the market variance left after Low risk\nthat is explained by the other 12 themes",
     "Share of\nexplained\nvariance", "output/figure1b_theme_importance_after_low_risk.png",
     paste0("Market and themes residualized on the Low risk theme, Feb 1991-Dec 2025.\nShare of the remaining variance explained by the 12 themes: ",
            paste(sprintf("%s %.0f%%", r2b$country, 100 * r2b$r2), collapse = ", "), "."))

cat("Partial R2 after Low risk:\n"); print(r2b)
cat("\nTop 3 themes after controlling for Low risk:\n")
print(partial[order(-share), head(.SD, 3), by = country][, .(country, theme, share = round(share, 2))])
