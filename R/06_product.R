# 9-11-bosqichlar. Mahsulot dizayni, tariflash va bazis riski
# (dissertatsiya 2.4 va 3.1-3.2 bo'limlari).

attach_province_index <- function(panel, province) {
  # Viloyat indeksi (vanom) - HOSILDORLIK (s/ga) anomaliyasi, tuman indeksi bilan bir xil o'lchovda.
  res <- detrend_series(province$year, province$yield_ts_ga)
  vy  <- data.frame(year = province$year, vanom = res$anomaly_pct)
  panel$vanom <- NULL
  merge(panel, vy, by = "year", sort = FALSE)
}

payout_area_yield <- function(anomaly_pct, deductible = DEDUCTIBLE) {
  # (3.2)-formula: to'lov = max(0, -anomaliya - ozod)
  pmax(0, -anomaly_pct - deductible)
}

pop_sd <- function(x) sqrt(mean((x - mean(x))^2))   # numpy std (ddof = 0) bilan bir xil

hedging_effectiveness <- function(loss, payout) {
  # (3.1)-formula: HE = [1 - sd(L - I) / sd(L)] * 100
  if (pop_sd(loss) == 0) return(NA_real_)
  (1 - pop_sd(loss - payout) / pop_sd(loss)) * 100
}

simulate_index <- function(panel, deductible = DEDUCTIBLE, index_level = c("district", "province")) {
  index_level <- match.arg(index_level)
  base <- if (index_level == "district") panel$anom else panel$vanom
  loss <- pmax(0, -panel$anom)
  pay  <- payout_area_yield(base, deductible)
  data.frame(district = panel$district, year = panel$year, loss = loss, payout = pay,
             net = loss - pay, stringsAsFactors = FALSE)
}

summarise_sim <- function(sim) {
  list(sof_tarif_pct        = mean(sim$payout),
       tolov_chastotasi_pct = mean(sim$payout > 0) * 100,
       HE_pct               = hedging_effectiveness(sim$loss, sim$payout),
       korrelyatsiya        = cor(sim$payout, sim$loss))
}

basis_risk_table <- function(panel, threshold = LOSS_THRESHOLD) {
  rows <- lapply(split(panel, panel$district), function(g) {
    actual <- g$anom < -threshold
    paid   <- g$vanom < -threshold
    data.frame(tuman = g$district[1], zarar_yillari = sum(actual), qoplangan = sum(actual & paid),
               qoplanmagan = sum(actual & !paid), asossiz_tolov = sum(!actual & paid))
  })
  df  <- do.call(rbind, rows)
  tot <- data.frame(tuman = "JAMI", t(colSums(df[, -1])))
  rbind(df, tot)
}

tariff <- function(pure_rate, loading = RISK_LOADING, subsidy = SUBSIDY_RATE, cap = SUBSIDY_CAP) {
  # Sof tarifdan fermer to'loviga. cap = NA - subsidiya chegarasiz (VM 206-son qarori).
  loaded <- pure_rate * (1 + loading)
  raw    <- subsidy * loaded
  sub    <- if (is.na(cap)) raw else min(raw, cap)
  list(sof = pure_rate, yuklamali = loaded, subsidiya = sub, fermer = loaded - sub,
       chegara_ishladi = !is.na(cap) && raw > cap)
}

high_risk_districts <- function(panel, sd_min = ZONE_SD_MIN, r_max = ZONE_R_MAX) {
  # 3.2.3: anomaliya SD > sd_min va viloyat anomaliyasi bilan r < r_max
  hit <- vapply(split(panel, panel$district), function(g) {
    sd(g$anom) > sd_min && cor(g$anom, g$vanom) < r_max
  }, logical(1))
  sort(names(hit)[hit])
}

zone_tariffs <- function(panel, deductible = DEDUCTIBLE, high_risk = HIGH_RISK_ZONE) {
  if (is.null(high_risk)) high_risk <- high_risk_districts(panel)
  sim   <- simulate_index(panel, deductible, "district")
  zones <- list("Barqaror zona"         = !sim$district %in% high_risk,
                "Yuqori xavf zonasi"    =  sim$district %in% high_risk,
                "Yagona viloyat tarifi" = rep(TRUE, nrow(sim)))
  do.call(rbind, lapply(names(zones), function(z) {
    g <- sim[zones[[z]], ]
    t <- tariff(mean(g$payout))
    data.frame(zona = z, tumanlar = length(unique(g$district)), sof = t$sof,
               yuklamali = t$yuklamali, subsidiya = t$subsidiya, fermer = t$fermer,
               HE_pct = hedging_effectiveness(g$loss, g$payout),
               chastota_pct = mean(g$payout > 0) * 100, stringsAsFactors = FALSE)
  }))
}
