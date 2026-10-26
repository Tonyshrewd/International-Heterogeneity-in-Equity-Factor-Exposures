# JKP Country-Month Panel — Data Codebook

FIN 864 project: *International Heterogeneity in Equity Factor Exposures*
Data built on 26 October 2026.

## Research team

- Anthony Lamini
- Pavel A. Solopov
- Maria Sultana
- Dr. Mohammad Najand (faculty supervisor)

## How to use

All data is included in this repository. The clean files in `data/` are ready to use. To load the main panel in R:

```r
panel <- readRDS("data/panel_vw_cap.rds")
```

To rebuild everything from the raw data (packages needed: `data.table`, `readxl`):

```
Rscript scripts/02_build_panel.R   # ~ 1 min: unzips raw/zips, cleans, checks, builds the panel
Rscript scripts/01_download.R      # optional: re-download the raw files from JKP (e.g. for a newer release)
```

Or open the script in RStudio, set the working directory to this folder, and click Source.

The raw JKP files are stored zipped in `raw/zips/` because each unzipped factor CSV (180–190 MB) is over GitHub's 100 MB file limit. The build script unzips them into `raw/csv/` automatically.

## Data license

The JKP Global Factor Data is licensed under CC BY-NC 4.0 (non-commercial use with attribution). Cite: Jensen, T. I., Kelly, B., and Pedersen, L. H. (2023). "Is There a Replication Crisis in Finance?" *Journal of Finance*, 78(5), 2465–2518.

## Source and release

| Item | Value |
| --- | --- |
| Source | JKP Global Factor Data, `https://jkpfactors-data.s3.amazonaws.com/public/` (the files behind jkpfactors.com/data) |
| Files | `[all_countries]_[mkt|all_themes|all_factors]_[monthly]_[vw_cap|vw|ew].zip` |
| CSV date stamp inside the zips | 14 April 2026 |
| Server upload date (Last-Modified) | 4 August 2026 |
| Last month of data | December 2025 |
| Release name | The JKP release notes give no version number or date; cite the three stamps above. |
| Download log | `raw/download_log.csv` (URL, size, server date for every file) |

## Folder layout

```
scripts/   01_download.R, 02_build_panel.R, theme_map.csv, country_meta.csv
raw/       zips/ (original JKP downloads, never edited), factor_details.xlsx, download_log.csv
           csv/ is created locally by the scripts and not stored on GitHub
data/      clean datasets (below)
output/    coverage tables and data checks
```

## Main dataset: `data/panel_vw_cap.csv` (also `.rds`)

One row per **country × month**. 36,395 rows, 95 countries. Same layout for `panel_vw` and `panel_ew` (robustness weightings).

| Column | Meaning |
| --- | --- |
| `country` | ISO-3 code (USA, JPN, DEU, GBR, AUS, …) |
| `country_name`, `msci` | Name and MSCI class (developed / emerging / frontier / standalone / not rated), from Table 10 of the JKP documentation |
| `date` | Last day of the month |
| `mkt` | Country market **excess** return, USD, capped value weighted. **Dependent variable.** |
| `mkt_n_stocks` | Number of stocks in the market portfolio that month |
| `accruals` … `value` | The 13 theme returns (long–short, USD). **Main regressors.** |
| `n_themes` | How many of the 13 themes are available that month |
| `complete` | TRUE if `mkt` and all 13 themes are available — the estimation sample for Model 1 |

**Units:** all returns are decimals per month (0.01 = 1%).

The 13 theme columns: `accruals`, `debt_issuance`, `investment`, `low_leverage`, `low_risk`, `momentum`, `profit_growth`, `profitability`, `quality`, `seasonality`, `short_term_reversal`, `size`, `value`.

## Other datasets

| File | Contents |
| --- | --- |
| `data/factors_long_<w>.rds` / `.csv.gz` | 153 individual factors, long format: `country, date, theme, factor, direction, n_stocks, n_stocks_min, ret`. For Model 3. |
| `data/themes_long.csv` | Theme returns, long format, all three weightings, with `n_factors` (factors averaged that month) |
| `data/market_long.csv` | Market returns, long format, all three weightings |
| `data/factor_dictionary.csv` | For each factor: theme, full name, original citation and sample, sign, whether the original paper found it significant |

## How JKP builds these numbers (from the documentation, verified where noted)

- **Factor:** within each country and month, stocks are sorted into thirds on the characteristic (breakpoints from non-micro stocks). Factor = top third minus bottom third. Each third uses capped value weights (market cap winsorized at the NYSE 80th percentile).
- **Sign:** each factor is long the side the original paper found to earn more. **Verified:** the `ret` column in the factor files is *already signed*. Do not multiply by `direction` again.
- **Theme:** the equal-weighted average of its factors that month. **Verified:** the average of the factor returns matches the published theme return to machine precision (max gap 2e-16) for every theme, country and month (`output/check_theme_reconstruction.csv`).
- **Minimum data:** a factor needs at least 5 stocks in each leg, and at least 60 months in a country to be published.
- **Theme map:** taken from Table 9 of the documentation (`scripts/theme_map.csv`). The `group` column in `factor_details.xlsx` is the Hou–Xue–Zhang category, *not* the JKP theme.

## Coverage summary (`output/coverage_by_country.csv`)

- 95 countries have market returns; 71 have theme and factor returns.
- **49 countries** have at least 120 complete months (market + all 13 themes): 23 developed, 22 emerging, 2 frontier, 2 other.
- Pilot countries (complete sample, no gaps):

| Country | Complete from | Complete months | Median stocks |
| --- | --- | --- | --- |
| USA | Nov 1951 | 890 | 4,245 |
| GBR | Apr 1986 | 477 | 1,589 |
| JPN | Apr 1986 | 477 | 3,609 |
| AUS | Nov 1990 | 422 | 1,418 |
| DEU | Feb 1991 | 419 | 728 |

- **Common-history window, 5 pilot countries:** Feb 1991 – Dec 2025, **419 months**, no gaps (`output/common_window_feasibility.csv`).

## Known issues to keep in mind

1. **Gaps.** 26 countries have missing months inside their market series (e.g., GRC, PRT, BRA, SAU, EGY). The column `longest_complete_run` in the coverage table shows the longest unbroken stretch; for rolling regressions use that, not `complete_months`.
2. **Thin markets.** Some eligible countries have very few stocks (e.g., IRL ~21, PRT ~43, COL ~50 in the last five years). With 5 stocks per leg as the minimum, their factor returns are noisy. Flag them in the sensitivity analysis.
3. **US history.** The US market starts in 1926 but all 13 themes only from Nov 1951, because the accounting-based factors need Compustat data.
4. **Currency.** Everything is in USD, so non-US returns include exchange-rate movements.
5. **Excess return.** The JKP site states returns are excess returns in USD. Sanity check: US market averages 0.73%/month with an annualized Sharpe ratio of 0.43 (1926–2025), consistent with an excess return.
