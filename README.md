# ARIMA Time Series Forecasting of Industrial Production Indices

This repository contains a project where **Industrial Production Indices (IPI)** for carpet and rug manufacturing in France were forecasted using **ARIMA models** in R. The project demonstrates time series analysis, model selection, forecasting, and visualization of trends.

---

## Project Overview

- **Objective:** Forecast monthly industrial production indices (NAF rev. 2, class 13.93) from January 2013 to March 2024.
- **Data:** Seasonally and calendar-adjusted indices (CVS-CJO) from INSEE, normalized to a base of 100 in 2021.
- **Tools:** R, `stats` package for ARIMA modeling, `ggplot2` for visualization.
- **Methodology:** Box-Jenkins framework for ARIMA modeling.

---

## Methodology

### 1. Data Exploration
- Plotted raw series to visualize trends and identify non-stationarity.
- Checked for linear trends and intercepts using regression on time.

### 2. Stationarity Testing
- Conducted **Augmented Dickey-Fuller (ADF)** and **KPSS tests**.
- Results indicated the series was non-stationary at level.
- Differentiated the series once to achieve stationarity (confirmed by ADF and KPSS tests on differenced series).

### 3. ARMA/ARIMA Modeling
- Examined **ACF** and **PACF** plots to select maximum orders `p` and `q`.
- Tested ARMA(p,q) models and validated using:
  - Parameter significance
  - Ljung-Box test for residual autocorrelation
  - AIC and BIC for model selection
- Selected the best-fit model: **ARIMA(1,1,1)**
  - AR(1) coefficient: 0.5016
  - MA(1) coefficient: -0.8196

### 4. Forecasting
- Produced forecasts for future indices (`X_{T+1}` and `X_{T+2}`).
- Computed confidence intervals based on Gaussian white-noise residual assumption.
- Visualized univariate forecasts and bivariate confidence ellipses.
- Observed high uncertainty due to residual variability and external shocks (e.g., COVID-19).

### 5. Open Question
- Investigated how additional real-time series `Y_t` could improve forecasts of `X_{T+1}` using Granger causality principles.

---

## Visualizations

1. Raw and differenced series  
2. ACF and PACF plots  
3. Forecasts with univariate and bivariate confidence regions  

*(Refer to the figures in the PDF or R scripts for full visualizations.)*


