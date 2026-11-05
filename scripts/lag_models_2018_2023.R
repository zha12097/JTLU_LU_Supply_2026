##############################################################################
#### Lag models 2018-2023: 1-year, 3-year and 5-year lag                  ####
#### Same specification as the main model in all three models;           ####
#### Change of Sale Price (SalePrice_CHG) is the only random parameter    ####
#### (normal, 200 Halton draws, panel). The three models differ only in   ####
#### the lag variable:                                                    ####
####   1-year lag: C_1Y   |   3-year lag: C_3Y   |   5-year lag: C_5Y      ####
####                                                                      ####
#### Produces Table 5 (2018-2023 columns) and Tables A.4-A.6.             ####
####                                                                      ####
#### Input : data/model_input_1y_2018_2023.csv                            ####
####         data/model_input_3y_2018_2023.csv                            ####
####         data/model_input_5y_2018_2023.csv                            ####
#### Output: outputs/lag_models_2018_2023/                                ####
####   Summary_<lag>.txt          full summary(model) report              ####
####   Coefficients_<lag>.csv     all coefficients                        ####
####   Model_<lag>.rds            fitted model                            ####
####   Lag_Comparison_Market_Signals.csv   market signals side by side     ####
####   Lag_Model_Fit.csv          log-likelihood, AIC, BIC, sample sizes   ####
#### Run from the repository root.                                        ####
##############################################################################

OUT_DIR <- file.path("outputs", "lag_models_2018_2023")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
op <- function(f) file.path(OUT_DIR, f)
R_DRAWS <- 200

library("mlogit")
library("tidyverse")

lags <- tribble(
  ~lag,          ~file,                                              ~lag_var,
  "1-Year Lag",  file.path("data", "model_input_1y_2018_2023.csv"),  "C_1Y",
  "3-Year Lag",  file.path("data", "model_input_3y_2018_2023.csv"),  "C_3Y",
  "5-Year Lag",  file.path("data", "model_input_5y_2018_2023.csv"),  "C_5Y")

# Model specification; the first term is the lag variable
make_formula <- function(lag_var) as.formula(paste(
  "choice ~", lag_var, "+ Rent_Adj + SalePrice_CHG + Cap_Rate_CHG + Lease_Deal_CHG + Sale_List |",
  "ParcelArea + CoStar_Entropy + a_AM_OFC + a_MD_IND + a_PM_RTL + EPOI_1km +",
  "ResAll_C_1km + BF_A_1km + BUID_C_100m + MJLC_RDS_L_1km + RAILRTS_A_1km + RAILRTS_L_100m +",
  "AIRP_DIST + BSTP_DIST + LOS_DIST + LU_COM + LU_IND + LU_RTL + LU_RES + Toronto_DTMid +",
  "Air + Bramp + CaleMilton + CoreDurham + Missi + OakBurling + Out_York + Outlying_Durham +",
  "Toronto_East + Toronto_West + York_South + OS_A_1km + NBSTP_DIST | 0"))

stars <- function(p) ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.10, "*", "")))
tag   <- function(lag) gsub("[^A-Za-z0-9]+", "", lag)          # "1-Year Lag" -> "1YearLag"

coefs <- list(); fits <- list()
for (i in seq_len(nrow(lags))) {
  L <- lags[i, ]
  cat("\n================ ", L$lag, " (", L$file, ") ================\n", sep = "")

  #### Data (same preparation as the main model) ####
  d <- read_csv(L$file, show_col_types = FALSE)
  d$id       <- as.integer(factor(d$GTHA_ID))
  d$choiceid <- as.integer(factor(10000 * d$GTHA_ID + d$Year))
  d$choice   <- as.logical(d$choice)
  d$Toronto_DTMid <- d$DT + d$Toronto_MidNorth
  d$BF_A_1km      <- d$BF_A_1km / 1E4
  d$EPOI_1km      <- d$EPOI_1km / 1E3
  d$ResAll_C_1km  <- d$ResAll_C_1km / 1E4
  d_idx <- dfidx(d, choice = "choice", idx = list(c("choiceid", "id"), "Dev_Type"))

  #### Estimation ####
  t0 <- Sys.time()
  m <- mlogit(make_formula(L$lag_var), data = d_idx,
              rpar = c(SalePrice_CHG = "n"), R = R_DRAWS, halton = NA, panel = TRUE)
  cat("Time taken:", round(difftime(Sys.time(), t0, units = "mins"), 1), "minutes\n")

  #### summary(model) report: console and text file ####
  s <- summary(m)
  print(s)
  capture.output(print(s), file = op(paste0("Summary_", tag(L$lag), ".txt")))
  saveRDS(m, op(paste0("Model_", tag(L$lag), ".rds")))

  #### Coefficients ####
  est <- coef(m); se <- sqrt(diag(as.matrix(vcov(m)))); z <- est / se; p <- 2 * pnorm(-abs(z))
  ct <- tibble(lag = L$lag, term = names(est),
               estimate = ifelse(grepl("^sd\\.", names(est)), abs(est), est),
               std_error = se, z_value = z, p_value = p, stars = stars(p))
  write_csv(ct, op(paste0("Coefficients_", tag(L$lag), ".csv")))
  coefs[[i]] <- ct

  LL <- as.numeric(logLik(m))[1]; K <- length(est); N <- n_distinct(d$choiceid)
  fits[[i]] <- tibble(lag = L$lag, choice_occasions = N, parcels = n_distinct(d$id),
                      developments = sum(d$choice & d$Dev_Type != "A_N_O"),
                      LogLik = LL, K = K, AIC = -2 * LL + 2 * K, BIC = -2 * LL + K * log(N))
}

#### Market signals side by side (Table 5); the lag variable shares one row ####
market <- bind_rows(coefs) |>
  mutate(term = ifelse(term %in% c("C_1Y", "C_3Y", "C_5Y"), "# Projects Developed n Years Ago", term)) |>
  filter(term %in% c("# Projects Developed n Years Ago", "Rent_Adj", "SalePrice_CHG", "Cap_Rate_CHG",
                     "Lease_Deal_CHG", "Sale_List", "sd.SalePrice_CHG")) |>
  select(lag, term, estimate, std_error, z_value, stars) |>
  pivot_wider(names_from = lag, values_from = c(estimate, std_error, z_value, stars), names_vary = "slowest")
write_csv(market, op("Lag_Comparison_Market_Signals.csv"))
write_csv(bind_rows(fits), op("Lag_Model_Fit.csv"))

cat("\n================ MARKET SIGNALS ACROSS LAG MODELS ================\n")
print(as.data.frame(market), digits = 4)
cat("\n================ MODEL FIT ================\n")
print(as.data.frame(bind_rows(fits)), digits = 7)
cat("\nFiles written to:", normalizePath(OUT_DIR), "\n")
