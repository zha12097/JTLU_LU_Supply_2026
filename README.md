# Land Use Supply at the Parcel Level in the Metropolitan-Scale Area: Framework, Application, and Behavioural Insights

Code and example data for:

> Zhang, J., Kiko, M., & Miller, E. J. *Land Use Supply at the Parcel Level in the Metropolitan-Scale Area: Framework, Application, and Behavioural Insights.* Submitted to *Journal of Transport and Land Use*.

## Data (`data/`)

Most input data are fully restricted by the data provider (CoStar Group) or available only to certain user groups (see the manuscript for details). The files below are therefore **synthetic**: they have the same structure as the real model inputs (columns, alternatives and choice sets), but every value is randomly generated, so results will not match the paper.

| File | Description |
|---|---|
| `model_input_1y_2015_2023.csv` | Parcel-year choice data, 1-year lag, 2015–2023 (main model) |
| `model_input_1y_2018_2023.csv` | Parcel-year choice data, 1-year lag, 2018–2023 |
| `model_input_3y_2018_2023.csv` | Parcel-year choice data, 3-year lag, 2018–2023 |
| `model_input_5y_2018_2023.csv` | Parcel-year choice data, 5-year lag, 2018–2023 |
| `distance_decay/model_input_*.csv` | Main-model data with the rent variable downscaled using alternative distance-decay functions (linear; exponential with half-life 1, 5, 10, 25 and 50 km) |

## Scripts (`scripts/`)

| Script | Description |
|---|---|
| `main_model_and_simulation.R` | Panel mixed logit (1-year lag, 2015–2023), MNL benchmark, model fit, prediction-success table and Monte Carlo simulation of market shares |
| `figure9_cumulative_share.R` | Observed vs. predicted cumulative market shares with a previous-year benchmark |
| `holdout_validation.R` | Five-fold parcel-level hold-out test of out-of-sample performance |
| `lag_models_2018_2023.R` | Models with 1-, 3- and 5-year lags estimated over 2018–2023 |
| `distance_decay_sensitivity.R` | Sensitivity of the market-signal estimates to the distance-decay function |
| `vif_test.R` | Multicollinearity (VIF) test of the numeric explanatory variables |

## License

MIT

## Acknowledgement

Claude Opus 5.5 (Anthropic) was used to tidy and restructure the code in this repository. The models, analyses and results are the authors' own.
