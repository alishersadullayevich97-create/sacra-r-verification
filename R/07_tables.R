# Dissertatsiyaning qo'shimcha jadvallari:
#   2.2.5  - qat'iy effektli panel regressiya (klaster SE) va ta'sir kattaligi
#   2.3.1  - 15 yillik oyna (2010-2024) bo'yicha Gumbel
#   2.3    - 2018-yil yo'qotishining kvantili va qaytarilish davri
#   2.3.3  - EQN asosidagi viloyat tariflari va burn-rate
#   2.4.1, 2.4.3, 3.2.2 - tumanlar bo'yicha risk, HE va sof tarif
#   3.1.2  - issiqlik indeksi (35 C+ kunlar) mahsuloti
#   3.1.4  - shartli ozod sezgirligi

# ---------- 2.2.5: panel regressiya ----------
fe_heat <- function(df, rhs) {
  # Klaster-robust SE (tuman bo'yicha): HC1 = (n-1)/(n-k) va G/(G-1) tuzatishlari;
  # p-qiymat normal taqsimotdan (df = Inf). Bu statsmodels'ning
  # cov_type = "cluster" sozlamasi bilan bir xil.
  d  <- df[complete.cases(df[, c("anom", "heat35")]), ]
  m  <- lm(as.formula(paste("anom ~", rhs)), data = d)
  V  <- sandwich::vcovCL(m, cluster = ~district, type = "HC1")
  ct <- lmtest::coeftest(m, vcov. = V, df = Inf)
  c(beta = ct["heat35", 1], p = ct["heat35", 4], n = nrow(d))
}

panel_fe_table <- function(panel, drop_year = 2018) {
  specs <- c("Tuman effektli"                 = "heat35 + factor(district)",
             "Tuman va yil effektli"          = "heat35 + factor(district) + factor(year)",
             "To'liq model (+ yog'in, namlik)" =
               "heat35 + precip_sum + rh_min_mean + factor(district) + factor(year)")
  do.call(rbind, lapply(names(specs), function(nm) {
    a <- fe_heat(panel, specs[[nm]])
    b <- fe_heat(panel[panel$year != drop_year, ], specs[[nm]])
    data.frame(spetsifikatsiya = nm, beta_barcha = a[["beta"]], p_barcha = a[["p"]],
               beta_2018siz = b[["beta"]], p_2018siz = b[["p"]],
               n = paste0(a[["n"]], "/", b[["n"]]), stringsAsFactors = FALSE)
  }))
}

within_sd <- function(panel, col) {
  u <- panel[[col]] - ave(panel[[col]], panel$district)
  sd(u - ave(u, panel$year))
}

effect_size <- function(panel, beta) {
  # 2.2.8: ta'sir = beta x within SD; r = beta x sd_x / sd_y
  sx <- within_sd(panel, "heat35"); sy <- within_sd(panel, "anom")
  list(within_sd_kun = sx, tasir_pct = beta * sx, r = beta * sx / sy, R2 = (beta * sx / sy)^2)
}

# ---------- 2.3: EQN qo'shimchalari ----------
window_evt <- function(province, start = 2010) {
  w   <- province[province$year >= start, ]
  res <- detrend_series(w$year, w$production_t)
  g   <- fit_gumbel(loss_series(res$anomaly_abs))
  mp  <- mean(w$production_t)
  list(g = g, rl = return_levels(g$mu, g$sigma, mean_production = mp),
       vc = var_cvar(g$mu, g$sigma, mean_production = mp))
}

event_return_period <- function(mu, sigma, loss, n_years) {
  q <- evd::pgumbel(loss, mu, sigma)
  T <- 1 / (1 - q)
  list(kvantil_pct = q * 100, davr_yil = T, P_n_yilda_pct = (1 - (1 - 1 / T)^n_years) * 100)
}

evt_tariffs <- function(mu, sigma, losses, mean_production, triggers = c(90, 92.5, 95),
                        loading = 0.25, subsidy = 0.5) {
  lp <- losses / mean_production * 100
  do.call(rbind, lapply(triggers, function(tr) {
    ded <- 100 - tr
    f   <- function(x) pmax(0, x / mean_production * 100 - ded) * evd::dgumbel(x, mu, sigma)
    ev  <- integrate(f, ded * mean_production / 100, evd::qgumbel(0.999999, mu, sigma),
                     subdivisions = 500L, rel.tol = 1e-10)$value
    data.frame(trigger_pct = tr, sof_EQN = ev, yuklamali = ev * (1 + loading),
               subsidiyadan_keyin = ev * (1 + loading) * (1 - subsidy),
               burn_rate = mean(pmax(0, lp - ded)))
  }))
}

# ---------- 2.4, 3.1, 3.2: tumanlar va mahsulot ----------
district_table <- function(panel, deductible = 10) {
  st <- simulate_index(panel, deductible, "district")
  sv <- simulate_index(panel, deductible, "province")
  out <- do.call(rbind, lapply(sort(unique(panel$district)), function(d) {
    g <- st[st$district == d, ]; gv <- sv[sv$district == d, ]; gp <- panel[panel$district == d, ]
    data.frame(tuman = d, n = nrow(gp), sd_pct = sd(gp$anom), r_viloyat = cor(gp$anom, gp$vanom),
               HE_tuman = hedging_effectiveness(g$loss, g$payout),
               HE_viloyat = hedging_effectiveness(gv$loss, gv$payout),
               sof_tarif = mean(g$payout), chastota_pct = mean(g$payout > 0) * 100,
               stringsAsFactors = FALSE)
  }))
  out[order(-out$sd_pct), ]
}

province_anomaly_sd <- function(province) {
  sd(detrend_series(province$year, province$yield_ts_ga)$anomaly_pct)
}

simulate_heat_index <- function(panel, beta, percentiles = c(50, 60, 70, 75, 80)) {
  # 3.1.2: to'lov = beta x max(0, issiq kunlar - tuman protsentili)
  # (protsentil: chiziqli interpolyatsiya, R type 7 = numpy "linear").
  loss <- pmax(0, -panel$anom)
  do.call(rbind, lapply(percentiles, function(q) {
    thr <- ave(panel$heat35, panel$district,
               FUN = function(s) quantile(s, q / 100, type = 7, names = FALSE))
    pay <- beta * pmax(0, panel$heat35 - thr)
    data.frame(protsentil = q, sof_tarif = mean(pay), chastota_pct = mean(pay > 0) * 100,
               HE_pct = hedging_effectiveness(loss, pay), korrelyatsiya = cor(pay, loss))
  }))
}

deductible_table <- function(panel, levels = c(5, 7.5, 10, 12.5, 15, 20)) {
  do.call(rbind, lapply(levels, function(d) {
    s <- summarise_sim(simulate_index(panel, d, "district"))
    data.frame(ozod_pct = d, sof_tarif_pct = s$sof_tarif_pct,
               tolov_chastotasi_pct = s$tolov_chastotasi_pct, HE_pct = s$HE_pct,
               korrelyatsiya = s$korrelyatsiya)
  }))
}
