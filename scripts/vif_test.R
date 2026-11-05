##############################################################################
#### Variance inflation factors for the numeric explanatory variables     ####
#### of the main model (no model estimation needed)                       ####
####                                                                      ####
#### Produces Table A.1.                                                  ####
####                                                                      ####
#### Input : data/model_input_1y_2015_2023.csv                            ####
#### Output: outputs/vif_test/VIF_Results.csv                             ####
#### Run from the repository root.                                        ####
##############################################################################

DATA_FILE <- file.path("data", "model_input_1y_2015_2023.csv")
OUT_DIR   <- file.path("outputs", "vif_test")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

##############################################
#### 1. DATA                              ####
##############################################
sample_data <- read.csv(DATA_FILE)

# Same scaling as in the main model
sample_data$BF_A_1km     <- sample_data$BF_A_1km / 1E4
sample_data$EPOI_1km     <- sample_data$EPOI_1km / 1E3
sample_data$ResAll_C_1km <- sample_data$ResAll_C_1km / 1E4

##############################################
#### 2. VARIABLES                         ####
##############################################
vif_predictors <- c(
  "C_1Y", "Rent_Adj", "SalePrice_CHG", "Cap_Rate_CHG", "Lease_Deal_CHG", "Sale_List",
  "ParcelArea", "CoStar_Entropy", "a_AM_OFC", "a_MD_IND", "a_PM_RTL",
  "EPOI_1km", "ResAll_C_1km", "BF_A_1km", "BUID_C_100m",
  "MJLC_RDS_L_1km", "RAILRTS_A_1km", "RAILRTS_L_100m",
  "AIRP_DIST", "BSTP_DIST", "LOS_DIST", "OS_A_1km", "NBSTP_DIST"
)

missing_predictors <- setdiff(vif_predictors, names(sample_data))
if (length(missing_predictors) > 0)
  stop("Variables missing from the dataset: ", paste(missing_predictors, collapse = ", "))

vif_data <- na.omit(as.data.frame(lapply(sample_data[, vif_predictors],
                                         function(x) as.numeric(as.character(x)))))
cat("Rows used:", nrow(vif_data), "| variables:", ncol(vif_data), "\n\n")

##############################################
#### 3. VIF                               ####
##############################################
# VIF_k = 1 / (1 - R2_k), where R2_k is from regressing variable k on all others
calculate_vif <- function(df) {
  vars <- names(df)
  results <- do.call(rbind, lapply(vars, function(target_var) {
    fit <- lm(reformulate(setdiff(vars, target_var), response = target_var), data = df)
    r2  <- summary(fit)$r.squared
    vif <- if (r2 >= 0.999999999) Inf else 1 / (1 - r2)
    data.frame(variable = target_var, r_squared = r2, vif = vif,
               tolerance = if (is.infinite(vif)) 0 else 1 / vif)
  }))
  results[order(-results$vif), ]
}

vif_results <- calculate_vif(vif_data)
print(vif_results, row.names = FALSE)

write.csv(vif_results, file.path(OUT_DIR, "VIF_Results.csv"), row.names = FALSE)
cat("\nVIF results written to:", normalizePath(OUT_DIR), "\n")
