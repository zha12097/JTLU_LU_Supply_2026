##############################################################################
#### Figure 9: observed vs. predicted cumulative market share by          ####
#### development type, GTHA 2015-2023 (1-year-lag model)                  ####
####                                                                      ####
#### Benchmark: previous-year persistence. The share of type j predicted  ####
#### for year t is the observed share of type j in year t-1, i.e. the     ####
#### same information set (outcomes up to t-1) as the one-step-ahead      ####
#### model simulation. 2015 has no previous year, so the benchmark starts ####
#### in 2016; its 2015 cumulative value is the observed 2015 share.       ####
####                                                                      ####
#### Input : outputs/main_model_and_simulation/Simulation_Shares.csv      ####
####         (written by main_model_and_simulation.R)                     ####
####         data/model_input_1y_2015_2023.csv                            ####
####         (only used to count choice occasions per year)              ####
#### Output: outputs/figure9_cumulative_share/                            ####
####         Figure9_Cumulative_Share.png/.pdf                            ####
####         Figure9_Benchmark_RMSE.csv  (model vs benchmark, 2016-2023)  ####
#### Run from the repository root.                                        ####
##############################################################################

library(dplyr)
library(ggplot2)

SIM_FILE  <- file.path("outputs", "main_model_and_simulation", "Simulation_Shares.csv")
DATA_FILE <- file.path("data", "model_input_1y_2015_2023.csv")
OUT_DIR   <- file.path("outputs", "figure9_cumulative_share")
EVAL_FROM <- 2016
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

sh <- read.csv(SIM_FILE)

# choice occasions per year (one chosen row per parcel-year)
occ <- read.csv(DATA_FILE)[, c("Year", "choice")] |>
  filter(choice == 1) |> count(Year, name = "n_occ")

#### 1. Previous-year persistence benchmark ####
yr <- sh |> filter(Measure == "Yearly share") |>
  select(Year, Alternative, observed) |>
  left_join(occ, by = "Year") |>
  arrange(Alternative, Year) |>
  group_by(Alternative) |>
  mutate(bench_yearly = lag(observed),
         # predicted counts: observed counts in 2015, persistence afterwards
         bench_count  = ifelse(is.na(bench_yearly), observed, bench_yearly) * n_occ,
         bench_cumul  = cumsum(bench_count) / cumsum(n_occ)) |>
  ungroup()

d_all <- sh |>
  left_join(yr |> select(Year, Alternative, bench_yearly, bench_cumul), by = c("Year", "Alternative")) |>
  mutate(benchmark = ifelse(Measure == "Yearly share", bench_yearly, bench_cumul)) |>
  select(-bench_yearly, -bench_cumul)

#### 2. RMSE (percentage points) and 95%-band coverage, 2016-2023 ####
ev <- d_all |> filter(Year >= EVAL_FROM)
rmse_tab <- bind_rows(
  ev |> group_by(Measure, Alternative),
  ev |> mutate(Alternative = "All types pooled") |> group_by(Measure, Alternative)) |>
  summarise(RMSE_model_pp     = sqrt(mean((mean - observed)^2)) * 100,
            RMSE_benchmark_pp = sqrt(mean((benchmark - observed)^2)) * 100,
            Coverage_95       = mean(observed >= lo & observed <= hi),
            .groups = "drop") |>
  mutate(RMSE_reduction_pct = (1 - RMSE_model_pp / RMSE_benchmark_pp) * 100) |>
  mutate(across(where(is.numeric), ~ round(.x, 3))) |>
  arrange(Measure, Alternative)
write.csv(rmse_tab, file.path(OUT_DIR, "Figure9_Benchmark_RMSE.csv"), row.names = FALSE)
print(as.data.frame(rmse_tab))

#### 3. Figure 9: cumulative market share ####
d <- d_all |> filter(Measure == "Cumulative share")

type_lab <- c(A_N_O = "No Development", Retail = "Retail", Industrial = "Industrial",
              Office = "Office", Mixed = "Mixed")
type_col <- c("No Development" = "#77AC1F", "Retail" = "#009E73", "Industrial" = "#0072B2",
              "Office" = "#D55E00", "Mixed" = "#CC79A7")
d$Type  <- factor(type_lab[d$Alternative], levels = type_lab)
d$Panel <- factor(ifelse(d$Alternative == "A_N_O", "Dominant Type (No Development)",
                         "Other Development Types"),
                  levels = c("Dominant Type (No Development)", "Other Development Types"))

src_lev <- c("Observed", "Predicted", "Benchmark (previous-year share)")
lines <- bind_rows(
  d |> transmute(Year, Type, Panel, value = observed, Source = "Observed"),
  d |> transmute(Year, Type, Panel, value = mean,     Source = "Predicted"),
  d |> filter(Year >= EVAL_FROM) |>
    transmute(Year, Type, Panel, value = benchmark, Source = "Benchmark (previous-year share)")) |>
  mutate(Source = factor(Source, levels = src_lev))
pts <- filter(lines, Source != "Benchmark (previous-year share)")

g <- ggplot() +
  geom_ribbon(data = d, aes(Year, ymin = lo, ymax = hi, fill = Type), alpha = 0.15, colour = NA) +
  geom_line(data = lines, aes(Year, value, colour = Type, linetype = Source), linewidth = 0.8) +
  geom_point(data = pts, aes(Year, value, colour = Type, shape = Source), size = 2, fill = "white") +
  scale_colour_manual(values = type_col, name = "Development Type") +
  scale_fill_manual(values = type_col, guide = "none") +
  scale_linetype_manual(values = c("Observed" = "solid", "Predicted" = "dashed",
                                   "Benchmark (previous-year share)" = "dotted"), name = "Data Source") +
  scale_shape_manual(values = c("Observed" = 16, "Predicted" = 24), guide = "none") +
  scale_x_continuous(breaks = 2015:2023) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~ Panel, ncol = 1, scales = "free_y") +
  labs(x = "Year", y = "Cumulative Market Share", title = "Cumulative Market Share",
       caption = paste("Shaded band: 95% range over 500 Monte Carlo runs.",
                       "Benchmark: observed share of the previous year (from 2016).")) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.background = element_rect(fill = "grey90", colour = "grey80"),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(fill = NA, colour = "grey80"),
        legend.position = "bottom", legend.box = "vertical",
        legend.title = element_text(face = "bold"),
        legend.key.width = unit(1.4, "cm")) +
  guides(linetype = guide_legend(order = 1, override.aes = list(colour = "black")),
         colour   = guide_legend(order = 2))

ggsave(file.path(OUT_DIR, "Figure9_Cumulative_Share.png"), g, width = 9, height = 7.5, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "Figure9_Cumulative_Share.pdf"), g, width = 9, height = 7.5)
