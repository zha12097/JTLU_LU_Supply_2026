##############################################################################
#### Main model: panel mixed logit, 1-year lag, 2015-2023                ####
####                                                                      ####
#### Change of Sale Price (SalePrice_CHG) is the only random parameter    ####
#### (normal, 200 Halton draws); all other coefficients are fixed.        ####
####                                                                      ####
#### 1. Data                                                              ####
#### 2. Estimation: mixed logit and MNL benchmark                         ####
#### 3. Model fit, LR test (mixed logit vs MNL), prediction-success table ####
#### 4. Coefficients with conventional and submarket-clustered SEs        ####
#### 5. Monte Carlo simulation of yearly and cumulative market shares     ####
####                                                                      ####
#### Produces Tables 3, 4, A.3 and Table 6 (Panel A), and the simulated   ####
#### shares used for Figure 9 (figure9_cumulative_share.R).               ####
####                                                                      ####
#### Input : data/model_input_1y_2015_2023.csv                            ####
#### Output: outputs/main_model_and_simulation/                           ####
#### Run from the repository root.                                        ####
##############################################################################

##############################################
#### 0. SETTINGS                          ####
##############################################
DATA_FILE   <- file.path("data", "model_input_1y_2015_2023.csv")
OUT_DIR     <- file.path("outputs", "main_model_and_simulation")
R_DRAWS     <- 200            # Halton draws for estimation
N_SIM       <- 500            # Monte Carlo runs for the simulation
RANDOM_VAR  <- "SalePrice_CHG"
NO_LABEL    <- "A_N_O"
CLUSTER_COL <- "MarketName"   # CoStar submarket (constant within parcel)

library("mlogit")
library("tidyverse")
set.seed(1)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
op  <- function(f) file.path(OUT_DIR, f)
say <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

##############################################
#### 1. DATA                              ####
##############################################
sample_data <- read_csv(DATA_FILE, show_col_types = FALSE)
sample_data$id       <- as.integer(factor(sample_data$GTHA_ID))
sample_data$choiceid <- as.integer(factor(10000 * sample_data$GTHA_ID + sample_data$Year))
sample_data$choice   <- as.logical(sample_data$choice)
sample_data$Toronto_DTMid <- sample_data$DT + sample_data$Toronto_MidNorth
sample_data$BF_A_1km      <- sample_data$BF_A_1km / 1E4
sample_data$EPOI_1km      <- sample_data$EPOI_1km / 1E3
sample_data$ResAll_C_1km  <- sample_data$ResAll_C_1km / 1E4
sample_data2 <- dfidx(sample_data, choice = 'choice', idx = list(c('choiceid', 'id'), 'Dev_Type'))

observed_choices <- sample_data |> filter(choice) |> distinct(choiceid, Year, Dev_Type)
cluster_of_occ   <- distinct(sample_data, choiceid, cluster = .data[[CLUSTER_COL]])
stopifnot(!any(duplicated(cluster_of_occ$choiceid)))

market_vars <- c("C_1Y", "Rent_Adj", "SalePrice_CHG", "Cap_Rate_CHG", "Lease_Deal_CHG", "Sale_List")
f <- choice ~ C_1Y + Rent_Adj + SalePrice_CHG + Cap_Rate_CHG +
  Lease_Deal_CHG + Sale_List |
  ParcelArea + CoStar_Entropy +
  a_AM_OFC + a_MD_IND + a_PM_RTL + EPOI_1km +
  ResAll_C_1km + BF_A_1km + BUID_C_100m +
  MJLC_RDS_L_1km + RAILRTS_A_1km + RAILRTS_L_100m +
  AIRP_DIST + BSTP_DIST + LOS_DIST + LU_COM +
  LU_IND + LU_RTL + LU_RES + Toronto_DTMid +
  Air + Bramp + CaleMilton + CoreDurham + Missi +
  OakBurling + Out_York + Outlying_Durham +
  Toronto_East + Toronto_West + York_South +
  OS_A_1km + NBSTP_DIST | 0

##############################################
#### HELPERS (compatible with mlogit 1.x and 2.x)
##############################################
mx_LL    <- function(m) if (inherits(m, "micsr")) as.numeric(logLik(m))[1] else as.numeric(m$logLik)[1]
mx_K     <- function(m) length(coef(m))
mx_vcov  <- function(m) {
  if (inherits(m, "micsr")) return(as.matrix(vcov(m)))
  V <- solve(-m$hessian); dimnames(V) <- list(names(coef(m)), names(coef(m))); V
}
mx_probs <- function(m) { P <- as.matrix(m$probabilities); P[is.na(P)] <- 0; P }
mx_occ   <- function(m) { o <- as.data.frame(idx(m))$choiceid; o[!duplicated(o)] }

# Cluster-robust covariance by submarket. mlogit stores one gradient row per choice
# occasion; for the panel model a parcel's rows add up to its score, and parcels
# never span two submarkets, so summing rows by submarket gives the cluster scores.
cluster_vcov <- function(m) {
  occ <- mx_occ(m); stopifnot(nrow(m$gradient) == length(occ))
  S <- rowsum(m$gradient, cluster_of_occ$cluster[match(occ, cluster_of_occ$choiceid)])
  G <- nrow(S); B <- mx_vcov(m)
  V <- B %*% (crossprod(S) * G / (G - 1)) %*% B
  dimnames(V) <- dimnames(B); attr(V, "G") <- G; V
}
stars <- function(p) ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.10, "*", "")))
coef_table <- function(m) {
  est <- coef(m); se <- sqrt(diag(mx_vcov(m)))
  Vc <- cluster_vcov(m); G <- attr(Vc, "G"); se_c <- sqrt(diag(Vc))
  tibble(term = names(est), estimate = ifelse(grepl("^sd\\.", names(est)), abs(est), est),
         std_error = se, z_value = est / se, p_value = 2 * pnorm(-abs(est / se)),
         stars = stars(p_value),
         se_cluster = se_c, t_cluster = est / se_c,
         p_cluster = 2 * pt(-abs(est / se_c), df = G - 1),     # t with G-1 df (few clusters)
         stars_cluster = stars(p_cluster), n_clusters = G)
}

##############################################
#### 2. ESTIMATION                        ####
##############################################
say("Estimating the MNL benchmark ...")
mnl <- mlogit(f, data = sample_data2)          # same model, no random parameter

say("Estimating the mixed logit (R = ", R_DRAWS, ") ...")
t0 <- Sys.time()
mxl <- mlogit(f, data = sample_data2, rpar = setNames("n", RANDOM_VAR),
              R = R_DRAWS, halton = NA, panel = TRUE)
say("Done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
saveRDS(mxl, op("MXL_model.rds"))
saveRDS(mnl, op("MNL_model.rds"))

sink(op("Results_Console_Log.txt"), split = TRUE)

##############################################
#### 3. FIT AND PREDICTION SUCCESS        ####
##############################################
n_occ <- n_distinct(sample_data$choiceid)
obs_counts <- table(observed_choices$Dev_Type)
LL_null <- sum(as.numeric(obs_counts) * log(as.numeric(obs_counts) / sum(obs_counts)))

LR <- 2 * (mx_LL(mxl) - mx_LL(mnl))
fit <- tibble(Model = c("Mixed logit", "MNL"),
              LogLik = c(mx_LL(mxl), mx_LL(mnl)), K = c(mx_K(mxl), mx_K(mnl))) |>
  mutate(Rho2 = 1 - LogLik / LL_null, Adj_Rho2 = 1 - (LogLik - K) / LL_null,
         LR_vs_null = 2 * (LogLik - LL_null),
         AIC = -2 * LogLik + 2 * K, BIC = -2 * LogLik + K * log(n_occ))
cat("\n=== MODEL FIT ===\n"); print(as.data.frame(fit), digits = 6)
cat("Null LL (aggregate shares):", round(LL_null, 3), "| choice occasions:", n_occ,
    "| parcels:", n_distinct(sample_data$id), "\n")
cat("LR test mixed logit vs MNL (H0: S.D. = 0): LR =", round(LR, 3),
    "| p (boundary-corrected) =", signif(0.5 * pchisq(LR, 1, lower.tail = FALSE), 4), "\n")

# Prediction-success table (McFadden 1977): rows = observed, columns = predicted probability mass
# PSI_j = N_jj / N_.j - N_j. / N_..  ; normalised PSI = PSI_j / (1 - N_j. / N_..)
P    <- mx_probs(mxl)
ids  <- mx_occ(mxl)
alts <- colnames(P)
obs  <- factor(observed_choices$Dev_Type[match(ids, observed_choices$choiceid)], levels = alts)
stopifnot(length(ids) == nrow(P), !any(is.na(obs)))
tab_raw <- rowsum(P, group = obs, reorder = FALSE)[alts, alts, drop = FALSE]
Ni <- rowSums(tab_raw); Nj <- colSums(tab_raw); N <- sum(tab_raw)
PSI <- diag(tab_raw) / Nj - Ni / N
pst <- data.frame(Observed = alts, round(100 * sweep(tab_raw, 1, Ni, "/"), 2),
                  Observed_N = round(Ni), Observed_Share = round(Ni / N, 4),
                  Predicted_Share = round(Nj / N, 4), PSI = round(PSI, 4),
                  PSI_Normalised = round(PSI / (1 - Ni / N), 4), check.names = FALSE)
cat("\n=== PREDICTION SUCCESS TABLE (%) ===\n"); print(pst)
write_csv(fit, op("Model_Fit.csv"))
write_csv(pst, op("Prediction_Success_Table.csv"))

##############################################
#### 4. COEFFICIENTS: CONVENTIONAL + CLUSTERED SEs
##############################################
coefs_mxl <- coef_table(mxl); coefs_mnl <- coef_table(mnl)
write_csv(coefs_mxl, op("Coefficients_MXL.csv"))
write_csv(coefs_mnl, op("Coefficients_MNL.csv"))
table3 <- coefs_mxl |> filter(term %in% c(market_vars, paste0("sd.", RANDOM_VAR)))
write_csv(table3, op("Table3_Market_Signals.csv"))
cat("\n=== MARKET SIGNALS: CONVENTIONAL vs CLUSTERED SEs ===\n")
print(as.data.frame(table3 |> select(term, estimate, std_error, p_value, stars,
                                     se_cluster, p_cluster, stars_cluster)), digits = 4)

##############################################
#### 5. SIMULATION OF MARKET SHARES       ####
#### - One-step-ahead: each year uses the observed lagged supply, market signals,
####   choice sets and parcel pool; simulated developments are not fed back.
#### - Probabilities are the model's fitted probabilities per choice occasion
####   (averaged over the distribution of the random coefficient).
#### - Each Monte Carlo run draws one alternative per occasion; N_SIM runs give
####   the mean and a 95% band. Occasions are linked to years by choiceid.
##############################################
cat("\n=== SIMULATION:", N_SIM, "Monte Carlo runs ===\n")
Pm  <- P / rowSums(P)
cpm <- t(apply(Pm, 1, cumsum)); cpm[, ncol(cpm)] <- 1
occ_year <- observed_choices$Year[match(ids, observed_choices$choiceid)]
stopifnot(!any(is.na(occ_year)))
years <- sort(unique(occ_year))

sim_counts <- array(0, dim = c(length(years), length(alts), N_SIM), dimnames = list(years, alts, NULL))
for (m in seq_len(N_SIM)) {
  pick <- max.col(cpm >= runif(nrow(cpm)), ties.method = "first")
  sim_counts[, , m] <- table(factor(occ_year, levels = years), factor(alts[pick], levels = alts))
}
check <- tibble(Alternative = alts, fitted = colMeans(Pm),
                simulated = apply(sim_counts, 2, sum) / N_SIM / nrow(Pm))
cat("\nCheck (overall predicted shares, should agree to about 3 decimals):\n"); print(check)
if (max(abs(check$fitted - check$simulated)) > 0.005) stop("Simulation does not reproduce the fitted shares.")

obs_counts_year <- table(factor(observed_choices$Year, levels = years),
                         factor(observed_choices$Dev_Type, levels = alts))
occ_per_year <- rowSums(obs_counts_year)
to_long <- function(cnt, type) {
  yearly <- sweep(cnt, 1, occ_per_year, "/")
  cumul  <- sweep(apply(cnt, 2, cumsum), 1, cumsum(occ_per_year), "/")
  bind_rows(as_tibble(yearly, rownames = "Year") |> pivot_longer(-Year, names_to = "Alternative") |> mutate(Measure = "Yearly share"),
            as_tibble(cumul,  rownames = "Year") |> pivot_longer(-Year, names_to = "Alternative") |> mutate(Measure = "Cumulative share")) |>
    mutate(Year = as.integer(Year), Series = type)
}
obs_long <- to_long(unclass(obs_counts_year), "Observed")
sim_long <- map_dfr(seq_len(N_SIM), function(m) to_long(sim_counts[, , m], "Model") |> mutate(run = m)) |>
  group_by(Year, Alternative, Measure) |>
  summarise(mean = mean(value), lo = quantile(value, 0.025), hi = quantile(value, 0.975), .groups = "drop")

sim_table <- obs_long |> select(Year, Alternative, Measure, observed = value) |>
  left_join(sim_long, by = c("Year", "Alternative", "Measure"))
write_csv(sim_table, op("Simulation_Shares.csv"))

cat("\nFiles written to:", normalizePath(OUT_DIR), "\n")
sink()
