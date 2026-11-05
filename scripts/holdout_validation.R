##############################################################################
#### Parcel-level hold-out test (5-fold), panel mixed logit 2015-2023     ####
#### Same specification as main_model_and_simulation.R                   ####
#### (SalePrice_CHG is the only random parameter).                        ####
####                                                                      ####
#### Protocol:                                                            ####
#### - Parcels are split into 5 folds. Each parcel's full 2015-2023       ####
####   sequence stays in one fold. Folds are balanced by each parcel's    ####
####   observed outcome (never developed / first development type).      ####
#### - The model is re-estimated on 4 folds and evaluated on the held-out ####
####   fold at the estimated coefficients; repeated for all 5 folds, so   ####
####   every parcel is predicted once by a model that did not use it.     ####
#### - Out-of-sample results are pooled over the 5 held-out folds and     ####
####   compared with the in-sample results of the full model on the same  ####
####   parcels.                                                           ####
#### - A parcel that is the only development of its type with a given     ####
####   0/1 dummy (e.g. the single Mixed case in Outlying Durham) is never ####
####   held out, because the coefficient could not be estimated without   ####
####   it; such parcels are listed in the console and not evaluated.      ####
####                                                                      ####
#### Produces Table 6 (Panel B).                                          ####
####                                                                      ####
#### Input : data/model_input_1y_2015_2023.csv                            ####
####         outputs/main_model_and_simulation/MXL_model.rds (full model, ####
####         from main_model_and_simulation.R; re-estimated if missing)   ####
#### Output: outputs/holdout_validation/                                  ####
####   HoldOut_Table.csv                     in-sample vs out-of-sample   ####
####   HoldOut_Fold_Summary.csv              per-fold check               ####
####   HoldOut_Prediction_Success_Table.csv  out-of-sample matrix         ####
#### Each fold is saved when finished; re-running skips finished folds.   ####
#### Run from the repository root.                                        ####
##############################################################################

##############################################
#### 0. SETTINGS                          ####
##############################################
DATA_FILE  <- file.path("data", "model_input_1y_2015_2023.csv")
OUT_DIR    <- file.path("outputs", "holdout_validation")
FULL_MODEL <- file.path("outputs", "main_model_and_simulation", "MXL_model.rds")
K_FOLDS    <- 5
R_DRAWS    <- 200        # same as the main model
RANDOM_VAR <- "SalePrice_CHG"
SEED       <- 2026

library("mlogit")
library("tidyverse")
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

make_dfidx <- function(d) dfidx(d, choice = 'choice', idx = list(c('choiceid', 'id'), 'Dev_Type'))
observed_choices <- sample_data |> filter(choice) |> distinct(choiceid, id, Year, Dev_Type)

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
#### 2. FOLDS (by parcel, balanced by outcome)
##############################################
NO_LABEL <- "A_N_O"
parcel_outcome <- observed_choices |>
  arrange(id, Year) |>
  group_by(id) |>
  summarise(outcome = { dev <- Dev_Type[Dev_Type != NO_LABEL]
  if (length(dev) == 0) "Never developed" else dev[1] }, .groups = "drop")
# Identification guard. Every alternative-specific coefficient on a 0/1 dummy
# (region, land use) needs at least one development of that type with the
# dummy = 1 in each training set, otherwise it is not estimable. So:
# (a) a parcel that is the ONLY such case (e.g. the single Mixed development
#     in Outlying Durham) is never held out (always in training);
# (b) folds are redrawn until no other such cell has all its parcels in one fold.
gen_vars <- all.vars(formula(Formula::Formula(f), lhs = 0, rhs = 2))
dummies  <- gen_vars[sapply(gen_vars, function(v) all(sample_data[[v]] %in% c(0, 1)))]
dev_rows <- sample_data |> filter(choice, Dev_Type != NO_LABEL)
cells <- map_dfr(dummies, function(v)
  dev_rows |> filter(.data[[v]] == 1) |> distinct(id, Dev_Type) |> mutate(dummy = v)) |>
  group_by(dummy, Dev_Type) |> mutate(n_parcels = n_distinct(id)) |> ungroup()
pinned <- unique(cells$id[cells$n_parcels == 1])

set.seed(SEED)
for (attempt in 1:1000) {
  folds <- parcel_outcome |>
    mutate(fold = 0L) |>                                   # 0 = never held out
    group_by(outcome) |>
    mutate(fold = ifelse(id %in% pinned, 0L,
                         sample(rep_len(sample(seq_len(K_FOLDS)), n())))) |>
    ungroup()
  bad <- cells |> filter(n_parcels > 1) |> left_join(folds, by = "id") |>
    group_by(dummy, Dev_Type) |>
    summarise(ok = n_distinct(fold[fold > 0]) > 1 | any(fold == 0), .groups = "drop") |>
    filter(!ok)
  if (nrow(bad) == 0) break
}
if (nrow(bad) > 0) stop("Could not find identifiable folds")
write_csv(folds, op("HoldOut_Fold_Assignment.csv"))
say("Fold assignment found on attempt ", attempt, ". Parcels by outcome and fold ",
    "(fold 0 = always in training, not evaluated):")
print(table(folds$outcome, folds$fold))
if (length(pinned) > 0) {
  cat("Parcels kept in every training set (only case of a dummy x alternative cell):\n")
  print(as.data.frame(cells |> filter(id %in% pinned, n_parcels == 1) |>
                        left_join(distinct(sample_data, id, GTHA_ID), by = "id") |>
                        select(GTHA_ID, Dev_Type, dummy)))
}
eval_ids <- folds$id[folds$fold > 0]

##############################################
#### 3. HELPERS                           ####
##############################################
# Per-occasion probabilities from an mlogit object, with choiceid and observed choice
get_probs <- function(m) {
  P <- as.matrix(fitted(m, type = "probabilities")); P[is.na(P)] <- 0
  ids <- as.data.frame(idx(m))$choiceid; ids <- ids[!duplicated(ids)]
  stopifnot(length(ids) == nrow(P))
  list(P = P, choiceid = ids,
       obs = observed_choices$Dev_Type[match(ids, observed_choices$choiceid)])
}

# Prediction-success table (McFadden 1977)
ps_table <- function(P, obs) {
  alts <- colnames(P); obs <- factor(obs, levels = alts)
  stopifnot(!any(is.na(obs)))
  tab <- rowsum(P, group = obs, reorder = FALSE)[alts, alts, drop = FALSE]
  Ni <- rowSums(tab); Nj <- colSums(tab); N <- sum(tab)
  PSI <- diag(tab) / Nj - Ni / N
  data.frame(Observed = alts, round(100 * sweep(tab, 1, Ni, "/"), 2),
             Observed_N = round(Ni), Observed_Share = round(Ni / N, 4),
             Predicted_Share = round(Nj / N, 4), PSI = round(PSI, 4),
             PSI_Normalised = round(PSI / (1 - Ni / N), 4), check.names = FALSE)
}

# Market-share null log-likelihood: shares from 'ref' applied to choices in 'obs'
LL_share_null <- function(obs, ref) {
  sh <- table(ref) / length(ref)
  sum(log(as.numeric(sh[as.character(obs)])))
}

code_msg <- c("1" = "converged (function values within tolerance)",
              "2" = "converged (gradient close to zero)",
              "3" = "converged (no higher value found; normal at the optimum)",
              "4" = "NOT CONVERGED: iteration limit reached")

##############################################
#### 4. IN-SAMPLE (full model)            ####
##############################################
if (file.exists(FULL_MODEL)) {
  say("Loading full model: ", FULL_MODEL)
  mxl_full <- readRDS(FULL_MODEL)
} else {
  say("Full model not found - estimating it on all parcels (R = ", R_DRAWS, ")")
  full_data <- make_dfidx(sample_data)
  mxl_full <- mlogit(f, data = full_data, rpar = setNames("n", RANDOM_VAR),
                     R = R_DRAWS, halton = NA, panel = TRUE)
  saveRDS(mxl_full, op("MXL_full_model.rds"))
}
# In-sample values on exactly the parcels that are evaluated out-of-sample
# (full-model coefficients, no re-estimation), so the comparison is like-for-like
all_obs  <- observed_choices$Dev_Type
mxl_eval <- update(mxl_full, start = coef(mxl_full),
                   data = make_dfidx(sample_data |> filter(id %in% eval_ids)),
                   iterlim = 0, print.level = 0)
ins <- get_probs(mxl_eval)
ins_LL      <- as.numeric(logLik(mxl_eval))
ins_LL_null <- LL_share_null(ins$obs, all_obs)
say("In-sample: LL = ", round(ins_LL, 2), " | occasions = ", length(ins$obs))

##############################################
#### 5. HOLD-OUT LOOP                     ####
##############################################
t_all <- Sys.time()
for (k in seq_len(K_FOLDS)) {
  fold_file <- op(sprintf("Fold_%d.rds", k))
  if (file.exists(fold_file)) { say("Fold ", k, " already done - skipped"); next }
  
  test_ids  <- folds$id[folds$fold == k]
  train_raw <- sample_data |> filter(!id %in% test_ids)
  test_raw  <- sample_data |> filter(id %in% test_ids)
  train_data <- make_dfidx(train_raw)
  test_data  <- make_dfidx(test_raw)
  say("Fold ", k, "/", K_FOLDS, ": estimating on ", n_distinct(train_raw$id), " parcels (",
      n_distinct(train_raw$choiceid), " occasions); holding out ", length(test_ids),
      " parcels (", n_distinct(test_raw$choiceid), " occasions)")
  
  t0 <- Sys.time()
  m_train <- mlogit(f, data = train_data, rpar = setNames("n", RANDOM_VAR),
                    R = R_DRAWS, halton = NA, panel = TRUE)
  conv <- as.character(m_train$est.stat$code)
  say("  estimated in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min; ",
      code_msg[conv], "; train LL = ", round(as.numeric(logLik(m_train)), 2))
  
  # Evaluate on the held-out parcels at the training coefficients (no re-estimation)
  m_test <- update(m_train, start = coef(m_train), data = test_data,
                   iterlim = 0, print.level = 0)
  stopifnot(isTRUE(all.equal(coef(m_test), coef(m_train))))
  te <- get_probs(m_test)
  tr <- get_probs(m_train)
  
  res <- list(fold = k, conv_code = conv,
              n_train_parcels = n_distinct(train_raw$id), n_test_parcels = length(test_ids),
              n_train_occ = length(tr$obs), n_test_occ = length(te$obs),
              LL_train = as.numeric(logLik(m_train)), LL_test = as.numeric(logLik(m_test)),
              LL_null_test = LL_share_null(te$obs, tr$obs),   # training shares applied to test
              coef = coef(m_train), test = te)
  saveRDS(res, fold_file)
  say("  held-out LL = ", round(res$LL_test, 2), " | LL per occasion: train ",
      round(res$LL_train / res$n_train_occ, 4), ", held-out ", round(res$LL_test / res$n_test_occ, 4))
  rm(m_train, m_test, train_data, test_data); gc(verbose = FALSE)
}
say("All folds finished. Total time: ", round(difftime(Sys.time(), t_all, units = "mins"), 1), " min")

##############################################
#### 6. RESULTS                           ####
##############################################
fold_res <- lapply(seq_len(K_FOLDS), function(k) readRDS(op(sprintf("Fold_%d.rds", k))))

fold_summary <- map_dfr(fold_res, function(r) tibble(
  Fold = r$fold, Convergence = unname(code_msg[r$conv_code]),
  Train_Parcels = r$n_train_parcels, Test_Parcels = r$n_test_parcels,
  Train_Occasions = r$n_train_occ, Test_Occasions = r$n_test_occ,
  Test_Developments = sum(r$test$obs != NO_LABEL),
  LL_per_Occ_Train = r$LL_train / r$n_train_occ,
  LL_per_Occ_Test  = r$LL_test / r$n_test_occ,
  Rho2_Test = 1 - r$LL_test / r$LL_null_test))
write_csv(fold_summary, op("HoldOut_Fold_Summary.csv"))
cat("\n=== PER-FOLD SUMMARY ===\n"); print(as.data.frame(fold_summary), digits = 4)
if (any(sapply(fold_res, `[[`, "conv_code") == "4"))
  warning("At least one fold did not converge - see HoldOut_Fold_Summary.csv")

# Pool the held-out predictions over folds
alts  <- colnames(ins$P)
P_oos <- do.call(rbind, lapply(fold_res, function(r) r$test$P[, alts, drop = FALSE]))
o_oos <- unlist(lapply(fold_res, function(r) r$test$obs))
stopifnot(length(o_oos) == length(ins$obs),
          setequal(unlist(lapply(fold_res, function(r) r$test$choiceid)), ins$choiceid))
oos_LL      <- sum(sapply(fold_res, `[[`, "LL_test"))
oos_LL_null <- sum(sapply(fold_res, `[[`, "LL_null_test"))

pst_in  <- ps_table(ins$P, ins$obs)
pst_oos <- ps_table(P_oos, o_oos)
write_csv(pst_oos, op("HoldOut_Prediction_Success_Table.csv"))
cat("\n=== OUT-OF-SAMPLE PREDICTION SUCCESS TABLE (%) ===\n"); print(pst_oos)

# The table for the paper: in-sample vs out-of-sample
lab <- c(A_N_O = "No development", Industrial = "Industrial", Mixed = "Mixed",
         Office = "Office", Retail = "Retail")
diag_in  <- diag(as.matrix(pst_in[, alts]));  diag_oos <- diag(as.matrix(pst_oos[, alts]))
holdout_table <- bind_rows(
  tibble(Measure = c("Choice occasions", "Log-likelihood", "Log-likelihood per occasion",
                     "Rho-squared (market-share null)"),
         In_sample      = c(length(ins$obs), ins_LL, ins_LL / length(ins$obs), 1 - ins_LL / ins_LL_null),
         Out_of_sample  = c(length(o_oos), oos_LL, oos_LL / length(o_oos), 1 - oos_LL / oos_LL_null)),
  tibble(Measure = paste0("Diagonal (%): ", lab[alts]),        In_sample = diag_in,  Out_of_sample = diag_oos),
  tibble(Measure = paste0("PSI: ", lab[alts]),                 In_sample = pst_in$PSI,  Out_of_sample = pst_oos$PSI),
  tibble(Measure = paste0("Normalised PSI: ", lab[alts]),      In_sample = pst_in$PSI_Normalised,
         Out_of_sample = pst_oos$PSI_Normalised)) |>
  mutate(Difference = Out_of_sample - In_sample)
write_csv(holdout_table, op("HoldOut_Table.csv"))
cat("\n=== HOLD-OUT TABLE: IN-SAMPLE vs OUT-OF-SAMPLE ===\n")
print(as.data.frame(holdout_table), digits = 4)

say("Done. Results in: ", normalizePath(OUT_DIR))