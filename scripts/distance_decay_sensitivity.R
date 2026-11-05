##############################################################################
#### Sensitivity of the market-signal coefficients to the distance-decay  ####
#### function used to downscale the rent variable                         ####
#### Same model in every column (SalePrice_CHG the only random parameter, ####
#### normal, 200 Halton draws);                                           ####
####   Squared (baseline) : model_input_1y_2015_2023.csv                  ####
####   Linear             : distance_decay/model_input_linear.csv         ####
####   Exponential, half-life h = 1, 5, 10, 25, 50 km:                     ####
####                        distance_decay/model_input_exp_hl<h>km.csv     ####
####                                                                      ####
####                                                                      ####
#### Produces Table A.2.                                                  ####
####                                                                      ####
#### Input : data/model_input_1y_2015_2023.csv                            ####
####         data/distance_decay/*.csv                                    ####
#### Output: outputs/distance_decay_sensitivity/                          ####
####         TableA2_Distance_Decay_Comparison.csv                        ####
#### Run from the repository root.                                        ####
##############################################################################

library("mlogit")
library("tidyverse")

DATA_DIR  <- "data"
DECAY_DIR <- file.path(DATA_DIR, "distance_decay")
OUT_DIR   <- file.path("outputs", "distance_decay_sensitivity")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

kernels <- tribble(
  ~label,      ~file,                                                ~rent_var,
  "Squared",   file.path(DATA_DIR,  "model_input_1y_2015_2023.csv"), "Rent_Adj",
  "Linear",    file.path(DECAY_DIR, "model_input_linear.csv"),       "Rent_LIN",
  "h = 1 km",  file.path(DECAY_DIR, "model_input_exp_hl01km.csv"),   "Rent_EXP",
  "h = 5 km",  file.path(DECAY_DIR, "model_input_exp_hl05km.csv"),   "Rent_EXP",
  "h = 10 km", file.path(DECAY_DIR, "model_input_exp_hl10km.csv"),   "Rent_EXP",
  "h = 25 km", file.path(DECAY_DIR, "model_input_exp_hl25km.csv"),   "Rent_EXP",
  "h = 50 km", file.path(DECAY_DIR, "model_input_exp_hl50km.csv"),   "Rent_EXP")

# Main model; "Rent" is replaced by each kernel's rent variable
f <- choice ~ C_1Y + Rent + SalePrice_CHG + Cap_Rate_CHG +
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

vars   <- c("C_1Y", "Rent", "SalePrice_CHG", "Cap_Rate_CHG", "Lease_Deal_CHG", "Sale_List")
labels <- c("# Projects Developed 1 Year Ago", "Rent Price", "Change of Sale Price",
            "Change of Cap Rate", "Change in # Leasing Deals", "# Sales Listed")
sig_stars <- function(p) ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.10, "*", "")))

results <- list()
for (i in seq_len(nrow(kernels))) {
  k <- kernels[i, ]
  message("Fitting ", k$label, " ...")
  d <- read_csv(k$file, show_col_types = FALSE)
  d$id       <- as.integer(factor(d$GTHA_ID))
  d$choiceid <- as.integer(factor(10000 * d$GTHA_ID + d$Year))
  d$choice   <- as.logical(d$choice)
  d$Toronto_DTMid <- d$DT + d$Toronto_MidNorth
  d$BF_A_1km      <- d$BF_A_1km / 1E4
  d$EPOI_1km      <- d$EPOI_1km / 1E3
  d$ResAll_C_1km  <- d$ResAll_C_1km / 1E4
  d$Rent <- d[[k$rent_var]]
  d <- d[!is.na(d$Rent), ]   # rows without a rent value: alternative unavailable
  d_idx <- dfidx(d, choice = "choice", idx = list(c("choiceid", "id"), "Dev_Type"))

  set.seed(1)
  m <- mlogit(f, data = d_idx, rpar = c(SalePrice_CHG = "n"), R = 200, halton = NA, panel = TRUE)

  est <- coef(m); se <- sqrt(diag(as.matrix(vcov(m)))); names(se) <- names(est)
  z <- est[vars] / se[vars]
  results[[k$label]] <- tibble(
    Variable = rep(labels, each = 4),
    Statistic = rep(c("Estimate", "Std. Error", "z-value", "Sig."), times = length(vars)),
    value = as.vector(rbind(sprintf("%.4f", est[vars]), sprintf("%.4f", se[vars]),
                            sprintf("%.4f", z), sig_stars(2 * pnorm(-abs(z))))))
}

table_A2 <- bind_rows(results, .id = "kernel") |>
  pivot_wider(names_from = kernel, values_from = value)
print(table_A2, n = Inf)
write_csv(table_A2, file.path(OUT_DIR, "TableA2_Distance_Decay_Comparison.csv"))
