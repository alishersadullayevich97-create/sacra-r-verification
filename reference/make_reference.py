"""Python (SACRA v1.1.0) natijalarini R bilan solishtirish uchun yagona jadvalga yozadi.

Ishga tushirish (SACRA repozitoriysining src/ papkasi ko'rsatiladi):

    python make_reference.py /yo'l/sacra-index-insurance/src

Natija: shu papkadagi python_v1.1.0.csv (ustunlar: key, value).
Har bir kalit R tomonda (run_all.R) aynan shu nom bilan hisoblanadi.
"""
import sys
from pathlib import Path

SRC = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(SRC))

import numpy as np
import pandas as pd
from scipy import stats

from config import DEDUCTIBLE
from s01_data import load_province, load_districts, load_climate, load_station
from s02_detrend import detrend_series, detrend_panel
from s03_climate import pearson_table, two_way_within
from s04_robustness import jackknife_correlation, permutation_test, partial_correlation
from s05_evt import (loss_series, fit_gumbel, check_gev_stability, gev_profile_likelihood,
                     return_levels, bootstrap_return_levels, var_cvar)
from s06_product import (attach_province_index, simulate_index, summarise, basis_risk_table,
                         zone_tariffs, high_risk_districts)
from s07_tables import (panel_fe_table, effect_size, window_evt, event_return_period,
                        evt_tariffs, district_table, province_anomaly_sd,
                        simulate_heat_index, deductible_table)

rows = []


def put(key, value):
    rows.append({"key": key, "value": value})


# ---------------- 1-3: ma'lumot va trend ----------------
v = load_province()
res = detrend_series(v.year, v.production_t)
v = v.assign(anom_abs=res["anomaly_abs"], anom=res["anomaly_pct"])
mp = v.production_t.mean()
i18 = list(v.year).index(2018)
put("province.n_years", len(v))
put("province.mean_t", mp)
put("province.cv_pct", v.production_t.std(ddof=1) / mp * 100)
put("province.trend_slope_t", res["slope"])
put("province.trend_p", res["p_value"])
put("province.anom2018_t", res["anomaly_abs"][i18])
put("province.anom2018_pct", res["anomaly_pct"][i18])

panel = detrend_panel(load_districts())
panel = attach_province_index(panel, v)
panel = panel.merge(load_climate(), on=["district", "year"], how="left")
put("panel.n_obs", len(panel))
put("panel.n_districts", panel.district.nunique())
put("panel.n_obs_no2018", int((panel.year != 2018).sum()))

# ---------------- 4: iqlim bog'liqligi ----------------
pt = pearson_table(panel, "anom", ["heat35", "heat30", "precip_sum", "rh_min_mean"])
for _, r in pt.iterrows():
    put(f"pearson.{r.omil}.r", r.r)
    put(f"pearson.{r.omil}.p", r.p)
w = two_way_within(panel, "anom", "heat35")
put("within.r", w["r"]); put("within.p", w["p"]); put("within.R2", w["r_squared"]); put("within.n", w["n"])

# ---------------- chang gipotezasi ----------------
st = load_station(); st = st[st.station == "Termiz"]
m = v[["year", "anom_abs"]].merge(st[["year", "dust_days", "heat_days"]], on="year")
m = m[(m.year >= 2014) & (m.year <= 2024)].reset_index(drop=True)
r0, p0, jk = jackknife_correlation(m.anom_abs, m.dust_days)
worst = jk.loc[jk.r.abs().idxmin()]
put("dust.n", len(m)); put("dust.r", r0); put("dust.p", p0)
put("dust.jackknife_worst_year", int(m.year[worst.chiqarilgan_indeks]))
put("dust.jackknife_worst_r", worst.r); put("dust.jackknife_worst_p", worst.p)
put("dust.spearman_rho", stats.spearmanr(m.anom_abs, m.dust_days)[0])
put("dust.permutation_p", permutation_test(m.anom_abs, m.dust_days))
pr, pp = partial_correlation(m.anom_abs, m.dust_days, m.heat_days)
put("dust.partial_r_heat", pr); put("dust.partial_p_heat", pp)
for _, r in jk.iterrows():
    put(f"dust.jackknife.{int(m.year[int(r.chiqarilgan_indeks)])}.r", r.r)

# ---------------- 5-8: EQN ----------------
L = loss_series(v.anom_abs)
gev = check_gev_stability(L)
put("gev.free_xi", gev["xi"]); put("gev.free_boot_ci_low", gev["ci_low"]); put("gev.free_boot_ci_high", gev["ci_high"])
pl = gev_profile_likelihood(L)
for k in ("xi", "ci_low", "ci_high", "ll_gev", "ll_gumbel", "LR", "p"):
    put(f"gev.profile.{k}", pl[k])
g = fit_gumbel(L)
put("gumbel.mu", g["mu"]); put("gumbel.sigma", g["sigma"])
for _, r in return_levels(g["dist"], mean_production=mp).iterrows():
    put(f"rl.T{int(r.qaytarilish_davri)}.t", r.yoqotish_t); put(f"rl.T{int(r.qaytarilish_davri)}.pct", r.yoqotish_pct)
for _, r in var_cvar(g["dist"], mean_production=mp).iterrows():
    a = int(round(r.ishonch * 100))
    put(f"var{a}.pct", r.VaR_pct); put(f"cvar{a}.pct", r.CVaR_pct)
for _, r in bootstrap_return_levels(L, mean_production=mp).iterrows():
    T = int(r.qaytarilish_davri)
    put(f"boot.T{T}.lo", r["quyi_2.5"]); put(f"boot.T{T}.hi", r["yuqori_97.5"])

# ---------------- 9-11: mahsulot ----------------
for lvl in ("district", "province"):
    s = summarise(simulate_index(panel, DEDUCTIBLE, lvl))
    for k, val in s.items():
        put(f"index.{lvl}.{k}", val)
br = basis_risk_table(panel); tot = br.iloc[-1]
put("basis.loss_events", tot.zarar_yillari); put("basis.covered", tot.qoplangan)
put("basis.uncovered", tot.qoplanmagan); put("basis.false_payouts", tot.asossiz_tolov)
hr = high_risk_districts(panel)
put("zone.high_risk", ";".join(hr))
for _, r in zone_tariffs(panel, high_risk=hr).iterrows():
    z = {"Barqaror zona": "stable", "Yuqori xavf zonasi": "high", "Yagona viloyat tarifi": "single"}[r.zona]
    for k in ("tumanlar", "sof", "yuklamali", "subsidiya", "fermer", "HE_pct", "chastota_pct"):
        put(f"zone.{z}.{k}", r[k])

# ---------------- s07: qo'shimcha jadvallar ----------------
fe = panel_fe_table(panel)
for i, tag in enumerate(("district_fe", "twoway_fe", "full_model")):
    for k in ("beta_barcha", "p_barcha", "beta_2018siz", "p_2018siz"):
        put(f"fe.{tag}.{k}", fe.loc[i, k])
b = abs(fe.loc[1, "beta_2018siz"])
for k, val in effect_size(panel, b).items():
    put(f"effect.{k}", val)
for k, val in event_return_period(g["dist"], L[i18], len(v)).items():
    put(f"event2018.{k}", val)
g15, rl15, vc15 = window_evt(v)
put("window2010.mu", g15["mu"]); put("window2010.sigma", g15["sigma"])
for _, r in rl15.iterrows():
    put(f"window2010.rl.T{int(r.qaytarilish_davri)}.pct", r.yoqotish_pct)
for _, r in vc15.iterrows():
    a = int(round(r.ishonch * 100)); put(f"window2010.cvar{a}.pct", r.CVaR_pct)
for _, r in evt_tariffs(g["dist"], L, mp).iterrows():
    t = str(r.trigger_pct)
    put(f"evt_tariff.{t}.sof_EQN", r.sof_EQN); put(f"evt_tariff.{t}.burn_rate", r.burn_rate)
put("province.yield_anom_sd", province_anomaly_sd(v))
for _, r in district_table(panel).iterrows():
    for k in ("n", "sd_pct", "r_viloyat", "HE_tuman", "HE_viloyat", "sof_tarif", "chastota_pct"):
        put(f"district.{r.tuman}.{k}", r[k])
for _, r in simulate_heat_index(panel, beta=b).iterrows():
    for k in ("sof_tarif", "chastota_pct", "HE_pct", "korrelyatsiya"):
        put(f"heat_index.P{int(r.protsentil)}.{k}", r[k])
for _, r in deductible_table(panel).iterrows():
    for k in ("sof_tarif_pct", "tolov_chastotasi_pct", "HE_pct", "korrelyatsiya"):
        put(f"deductible.{r.ozod_pct}.{k}", r[k])

out = Path(__file__).resolve().parent / "python_v1.1.0.csv"
pd.DataFrame(rows).to_csv(out, index=False)
print(f"{len(rows)} ta qiymat yozildi: {out}")
