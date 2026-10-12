# 2-3-bosqichlar. Texnologik trendni ajratish.
# Anomaliya (dissertatsiya (2.1)-formulasi): A_it = (Y_it - Yhat_it) / Yhat_it * 100

detrend_series <- function(years, values) {
  fit   <- lm(values ~ years)
  cf    <- summary(fit)$coefficients
  trend <- as.numeric(fitted(fit))
  list(slope       = unname(coef(fit)[2]),
       intercept   = unname(coef(fit)[1]),
       p_value     = cf[2, 4],
       r_squared   = summary(fit)$r.squared,
       trend       = trend,
       anomaly_pct = (values - trend) / trend * 100,
       anomaly_abs = values - trend)
}

detrend_panel <- function(df, value_col = "yield_ts_ga", group_col = "district") {
  # Har bir tuman uchun ALOHIDA trend (texnologik o'sish sur'ati tumanlarda har xil).
  parts <- lapply(split(df, df[[group_col]]), function(g) {
    g <- g[order(g$year), ]
    if (nrow(g) < MIN_YEARS_TREND) return(NULL)
    r <- detrend_series(g$year, g[[value_col]])
    g$trend <- r$trend
    g$anom  <- r$anomaly_pct
    g
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}
