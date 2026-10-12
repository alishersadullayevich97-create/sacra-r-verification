# SACRA v1.1.0 hisob-kitoblarining R'dagi mustaqil takrori.
#
# Ishga tushirish:
#   RStudio:  sacra-r.Rproj ni oching, so'ng  source("run_all.R")
#   Terminal: Rscript run_all.R
#
# Natijalar outputs/ papkasiga yoziladi. Python bilan solishtirish: source("compare.R")

find_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f) == 1) return(dirname(normalizePath(f)))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(dirname(normalizePath(of)))
  }
  if (file.exists(file.path("R", "00_config.R"))) return(normalizePath("."))
  stop("Repozitoriy ildizini topib bo'lmadi: sacra-r.Rproj ni oching yoki setwd() qiling")
}

ROOT <- find_root()
for (p in c("evd", "sandwich", "lmtest")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop(sprintf("'%s' paketi o'rnatilmagan: install.packages(c(\"evd\", \"sandwich\", \"lmtest\"))", p))
}
for (f in sort(list.files(file.path(ROOT, "R"), pattern = "\\.R$", full.names = TRUE))) source(f)

results <- list()
put <- function(key, value) results[[key]] <<- value
hdr <- function(t) cat("\n", strrep("=", 68), "\n", t, "\n", strrep("=", 68), "\n", sep = "")

# ---------------- 1-3: ma'lumot va trend ----------------
hdr("1-3. MA'LUMOT VA TREND")
v   <- load_province()
res <- detrend_series(v$year, v$production_t)
v$anom_abs <- res$anomaly_abs; v$anom <- res$anomaly_pct
mp  <- mean(v$production_t)
i18 <- which(v$year == 2018)
cat(sprintf("Viloyat: %d yil, o'rtacha %s t, CV %.1f%%\n", nrow(v),
            format(round(mp), big.mark = " "), sd(v$production_t) / mp * 100))
cat(sprintf("Trend: %+.0f t/yil (p = %.4f)\n", res$slope, res$p_value))
cat(sprintf("2018: %.0f t (%.1f%%)\n", res$anomaly_abs[i18], res$anomaly_pct[i18]))
put("province.n_years", nrow(v)); put("province.mean_t", mp)
put("province.cv_pct", sd(v$production_t) / mp * 100)
put("province.trend_slope_t", res$slope); put("province.trend_p", res$p_value)
put("province.anom2018_t", res$anomaly_abs[i18]); put("province.anom2018_pct", res$anomaly_pct[i18])

panel <- detrend_panel(load_districts())
panel <- attach_province_index(panel, v)
panel <- merge(panel, load_climate(), by = c("district", "year"), all.x = TRUE)
panel <- panel[order(panel$district, panel$year), ]; rownames(panel) <- NULL
cat(sprintf("Panel: %d kuzatuv, %d tuman\n", nrow(panel), length(unique(panel$district))))
put("panel.n_obs", nrow(panel)); put("panel.n_districts", length(unique(panel$district)))
put("panel.n_obs_no2018", sum(panel$year != 2018))

# ---------------- 4: iqlim bog'liqligi ----------------
hdr("4. IQLIM BOG'LIQLIGI")
pt <- pearson_table(panel, "anom", c("heat35", "heat30", "precip_sum", "rh_min_mean"))
print(transform(pt, r = round(r, 3), p = round(p, 4)), row.names = FALSE)
for (i in seq_len(nrow(pt))) { put(sprintf("pearson.%s.r", pt$omil[i]), pt$r[i]); put(sprintf("pearson.%s.p", pt$omil[i]), pt$p[i]) }
w <- two_way_within(panel, "anom", "heat35")
cat(sprintf("\nTuman+yil effektlaridan keyin: r = %+.3f, p = %.4f, R2 = %.3f\n", w$r, w$p, w$r_squared))
put("within.r", w$r); put("within.p", w$p); put("within.R2", w$r_squared); put("within.n", w$n)

# ---------------- chang gipotezasi ----------------
hdr("BARQARORLIK SINOVLARI (chang gipotezasi)")
st <- load_station(); st <- st[st$station == "Termiz", ]
m  <- merge(v[, c("year", "anom_abs")], st[, c("year", "dust_days", "heat_days")], by = "year")
m  <- m[m$year >= 2014 & m$year <= 2024, ]; m <- m[order(m$year), ]
jk <- jackknife_correlation(m$anom_abs, m$dust_days)
wi <- which.min(abs(jk$table$r))
rho  <- cor(m$anom_abs, m$dust_days, method = "spearman")
perm <- permutation_test(m$anom_abs, m$dust_days)
pc   <- partial_correlation(m$anom_abs, m$dust_days, m$heat_days)
cat(sprintf("Asosiy: r = %+.3f, p = %.4f, n = %d\n", jk$r0, jk$p0, nrow(m)))
cat(sprintf("Eng ta'sirchan kuzatuv: %d-yil -> r = %+.3f, p = %.4f\n", m$year[wi], jk$table$r[wi], jk$table$p[wi]))
cat(sprintf("Spearman: rho = %+.3f\nPermutatsiya p = %.4f\n", rho, perm))
cat(sprintf("Qisman (heat): r = %+.3f, p = %.4f (n-2), p = %.4f (n-3)\n", pc$r, pc$p_python, pc$p_dof))
put("dust.n", nrow(m)); put("dust.r", jk$r0); put("dust.p", jk$p0)
put("dust.jackknife_worst_year", m$year[wi]); put("dust.jackknife_worst_r", jk$table$r[wi])
put("dust.jackknife_worst_p", jk$table$p[wi]); put("dust.spearman_rho", rho); put("dust.permutation_p", perm)
put("dust.partial_r_heat", pc$r); put("dust.partial_p_heat", pc$p_python)
for (i in seq_len(nrow(m))) put(sprintf("dust.jackknife.%d.r", m$year[i]), jk$table$r[i])
write.csv(transform(jk$table, yil = m$year), file.path(OUT, "jackknife.csv"), row.names = FALSE)

# ---------------- 5-8: EQN ----------------
hdr("5-8. EKSTREMAL QIYMATLAR NAZARIYASI")
L   <- loss_series(v$anom_abs)
gev <- check_gev_stability(L)
cat(sprintf("Erkin GEV: xi = %+.3f, bootstrap CI = [%+.2f, %+.2f] (%d/%d muvaffaqiyatli)\n",
            gev$xi, gev$ci_low, gev$ci_high, gev$n_ok, gev$n_boot))
pl <- gev_profile_likelihood(L)
cat(sprintf("Profil ehtimollik: xi = %+.3f, 95%% CI = [%+.3f, %+.3f]; LR = %.3f, p = %.3f (log-L: GEV %.2f, Gumbel %.2f)\n",
            pl$xi, pl$ci_low, pl$ci_high, pl$LR, pl$p, pl$ll_gev, pl$ll_gumbel))
g <- fit_gumbel(L)
cat(sprintf("Gumbel: mu = %.0f, sigma = %.0f\n\n", g$mu, g$sigma))
rl <- return_levels(g$mu, g$sigma, mean_production = mp)
vc <- var_cvar(g$mu, g$sigma, mean_production = mp)
bs <- bootstrap_return_levels(L, mean_production = mp)
print(transform(rl, yoqotish_t = round(yoqotish_t), yoqotish_pct = round(yoqotish_pct, 1)), row.names = FALSE)
print(transform(vc, VaR_t = round(VaR_t), CVaR_t = round(CVaR_t), VaR_pct = round(VaR_pct, 1), CVaR_pct = round(CVaR_pct, 1)), row.names = FALSE)
print(data.frame(lapply(bs, function(x) round(x, 1)), check.names = FALSE), row.names = FALSE)
put("gev.free_xi", gev$xi); put("gev.free_boot_ci_low", gev$ci_low); put("gev.free_boot_ci_high", gev$ci_high)
for (k in c("xi", "ci_low", "ci_high", "ll_gev", "ll_gumbel", "LR", "p")) put(paste0("gev.profile.", k), pl[[k]])
put("gumbel.mu", g$mu); put("gumbel.sigma", g$sigma)
for (i in seq_len(nrow(rl))) {
  T <- rl$qaytarilish_davri[i]
  put(sprintf("rl.T%d.t", T), rl$yoqotish_t[i]); put(sprintf("rl.T%d.pct", T), rl$yoqotish_pct[i])
  put(sprintf("boot.T%d.lo", T), bs$quyi_2.5[i]); put(sprintf("boot.T%d.hi", T), bs$yuqori_97.5[i])
}
for (i in seq_len(nrow(vc))) {
  a <- round(vc$ishonch[i] * 100)
  put(sprintf("var%d.pct", a), vc$VaR_pct[i]); put(sprintf("cvar%d.pct", a), vc$CVaR_pct[i])
}
write.csv(rl, file.path(OUT, "return_levels.csv"), row.names = FALSE)
write.csv(bs, file.path(OUT, "bootstrap_ci.csv"), row.names = FALSE)
write.csv(as.data.frame(pl), file.path(OUT, "gev_profile.csv"), row.names = FALSE)

# ---------------- 9-11: mahsulot ----------------
hdr("9-11. MAHSULOT DIZAYNI VA BAZIS RISKI")
for (lvl in c("district", "province")) {
  s <- summarise_sim(simulate_index(panel, DEDUCTIBLE, lvl))
  cat(sprintf("%-9s indeksi: sof tarif %.2f%%, chastota %.0f%%, HE %.1f%%, r = %.2f\n",
              lvl, s$sof_tarif_pct, s$tolov_chastotasi_pct, s$HE_pct, s$korrelyatsiya))
  for (k in names(s)) put(sprintf("index.%s.%s", lvl, k), s[[k]])
}
br  <- basis_risk_table(panel); tot <- br[br$tuman == "JAMI", ]
cat(sprintf("\nBazis riski: %d ta zarar holati, %d tasi qoplangan (%.0f%%)\n",
            tot$zarar_yillari, tot$qoplangan, tot$qoplangan / tot$zarar_yillari * 100))
put("basis.loss_events", tot$zarar_yillari); put("basis.covered", tot$qoplangan)
put("basis.uncovered", tot$qoplanmagan); put("basis.false_payouts", tot$asossiz_tolov)
write.csv(br, file.path(OUT, "basis_risk.csv"), row.names = FALSE)

hr <- high_risk_districts(panel)
zt <- zone_tariffs(panel, high_risk = hr)
cat(sprintf("\nYuqori xavf zonasi (SD > 14%%, r < 0,5): %s\n", paste(hr, collapse = ", ")))
print(transform(zt, sof = round(sof, 2), yuklamali = round(yuklamali, 2), subsidiya = round(subsidiya, 2),
                fermer = round(fermer, 2), HE_pct = round(HE_pct, 1), chastota_pct = round(chastota_pct)),
      row.names = FALSE)
put("zone.high_risk", paste(hr, collapse = ";"))
zkey <- c("Barqaror zona" = "stable", "Yuqori xavf zonasi" = "high", "Yagona viloyat tarifi" = "single")
for (i in seq_len(nrow(zt))) for (k in c("tumanlar", "sof", "yuklamali", "subsidiya", "fermer", "HE_pct", "chastota_pct"))
  put(sprintf("zone.%s.%s", zkey[[zt$zona[i]]], k), zt[[k]][i])
write.csv(zt, file.path(OUT, "zone_tariffs.csv"), row.names = FALSE)

# ---------------- qo'shimcha jadvallar ----------------
hdr("QO'SHIMCHA JADVALLAR: 2.2.5, 2.3, 2.4, 3.1")
fe <- panel_fe_table(panel)
print(transform(fe, beta_barcha = round(beta_barcha, 3), p_barcha = round(p_barcha, 3),
                beta_2018siz = round(beta_2018siz, 3), p_2018siz = round(p_2018siz, 3)), row.names = FALSE)
ftag <- c("district_fe", "twoway_fe", "full_model")
for (i in 1:3) for (k in c("beta_barcha", "p_barcha", "beta_2018siz", "p_2018siz"))
  put(sprintf("fe.%s.%s", ftag[i], k), fe[[k]][i])
# Python 2.2.8 va 3.1.2 da 3 xonagacha yaxlitlangan koeffitsientni ishlatadi - shu takrorlanadi.
b  <- round(abs(fe$beta_2018siz[2]), 3)
es <- effect_size(panel, b)
cat(sprintf("\nTa'sir kattaligi (beta = %.3f): within SD = %.2f kun, ta'sir = %.2f%%, r = %.3f, R2 = %.3f\n",
            b, es$within_sd_kun, es$tasir_pct, es$r, es$R2))
for (k in names(es)) put(paste0("effect.", k), es[[k]])
ev <- event_return_period(g$mu, g$sigma, L[i18], nrow(v))
cat(sprintf("2018-yil: kvantil %.2f%%, qaytarilish davri ~%.0f yil, %d yilda kuzatilish ehtimoli %.1f%%\n",
            ev$kvantil_pct, ev$davr_yil, nrow(v), ev$P_n_yilda_pct))
for (k in names(ev)) put(paste0("event2018.", k), ev[[k]])

w15 <- window_evt(v)
cat(sprintf("\n2.3.1 - 15 yillik oyna (2010-2024): mu = %.0f, sigma = %.0f\n", w15$g$mu, w15$g$sigma))
put("window2010.mu", w15$g$mu); put("window2010.sigma", w15$g$sigma)
for (i in seq_len(nrow(w15$rl))) put(sprintf("window2010.rl.T%d.pct", w15$rl$qaytarilish_davri[i]), w15$rl$yoqotish_pct[i])
for (i in seq_len(nrow(w15$vc))) put(sprintf("window2010.cvar%d.pct", round(w15$vc$ishonch[i] * 100)), w15$vc$CVaR_pct[i])

et <- evt_tariffs(g$mu, g$sigma, L, mp)
print(data.frame(lapply(et, function(x) round(x, 2))), row.names = FALSE)
for (i in seq_len(nrow(et))) {
  tr <- sprintf("%.1f", et$trigger_pct[i])
  put(sprintf("evt_tariff.%s.sof_EQN", tr), et$sof_EQN[i]); put(sprintf("evt_tariff.%s.burn_rate", tr), et$burn_rate[i])
}
put("province.yield_anom_sd", province_anomaly_sd(v))

dt <- district_table(panel)
print(transform(dt, sd_pct = round(sd_pct, 1), r_viloyat = round(r_viloyat, 2), HE_tuman = round(HE_tuman, 1),
                HE_viloyat = round(HE_viloyat, 1), sof_tarif = round(sof_tarif, 2), chastota_pct = round(chastota_pct)),
      row.names = FALSE)
for (i in seq_len(nrow(dt))) for (k in c("n", "sd_pct", "r_viloyat", "HE_tuman", "HE_viloyat", "sof_tarif", "chastota_pct"))
  put(sprintf("district.%s.%s", dt$tuman[i], k), dt[[k]][i])

hi <- simulate_heat_index(panel, beta = b)
for (i in seq_len(nrow(hi))) for (k in c("sof_tarif", "chastota_pct", "HE_pct", "korrelyatsiya"))
  put(sprintf("heat_index.P%d.%s", hi$protsentil[i], k), hi[[k]][i])
dd <- deductible_table(panel)
for (i in seq_len(nrow(dd))) for (k in c("sof_tarif_pct", "tolov_chastotasi_pct", "HE_pct", "korrelyatsiya"))
  put(sprintf("deductible.%.1f.%s", dd$ozod_pct[i], k), dd[[k]][i])

for (nm in c("panel_fe", "districts", "evt_tariffs", "heat_index_product", "deductible")) {
  df <- switch(nm, panel_fe = fe, districts = dt, evt_tariffs = et, heat_index_product = hi, deductible = dd)
  write.csv(df, file.path(OUT, paste0(nm, ".csv")), row.names = FALSE)
}
write.csv(w15$rl, file.path(OUT, "evt_window_2010.csv"), row.names = FALSE)

kv <- data.frame(key = names(results),
                 value = vapply(results, function(x) if (is.character(x)) x else formatC(x, digits = 15, format = "g"),
                                character(1)), stringsAsFactors = FALSE)
write.csv(kv, file.path(OUT, "r_results.csv"), row.names = FALSE)
cat(sprintf("\n%d ta ko'rsatkich hisoblandi. Natijalar: %s\n", nrow(kv), OUT))
cat("Python bilan solishtirish: Rscript compare.R  (yoki source(\"compare.R\"))\n")
