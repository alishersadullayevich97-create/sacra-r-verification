# SACRA v1.1.0 — R'da mustaqil takrorlash

Dissertatsiyaning empirik hisob-kitoblari (SACRA v1.1.0, Python) R'da qaytadan
yozilgan va natijalar Python bilan avtomatik solishtiriladi. Maqsad — natijalar
dasturiy ta'minotga bog'liq emasligini ko'rsatish: boshqa til, boshqa kutubxonalar,
boshqa optimallashtirish algoritmlari bilan ham bir xil raqamlar chiqishi kerak.

- Asl tahlil: [sacra-index-insurance](https://github.com/alishersadullayevich97-create/sacra-index-insurance), teg `v1.1.0`
- Guvohnoma: DGU № 60955
- Muallif: A.S. Abdullayev

Bu repozitoriy yopiq (private) — tekshiruv uchun.

---

## Tez boshlash

R 4.2 yoki undan yangi versiya kerak. Paketlarni bir marta o'rnating:

```r
install.packages(c("evd", "sandwich", "lmtest"))
```

Sinovdan o'tgan muhit: R 4.3.3; evd 2.3.6.1; sandwich 3.1.0; lmtest 0.9.40.

**RStudio:** `sacra-r.Rproj` ni oching, so'ng konsolda:

```r
source("run_all.R")    # ~10 soniya; natijalar outputs/ papkasida
source("compare.R")    # Python bilan solishtirish jadvali
```

**Terminal:**

```
Rscript run_all.R
Rscript compare.R
```

## Tuzilish

| Fayl | Python moduli | Mazmuni |
|---|---|---|
| `R/00_config.R` | `config.py` | Parametrlar va ma'lumot sifati qoidalari |
| `R/01_data.R` | `s01_data.py` | Yuklash; 2006, Bandixon 2011–2019, 2010 o'tish yili, Qiziriq 2024 qoidalari |
| `R/02_detrend.R` | `s02_detrend.py` | Chiziqli trend, foizli anomaliya |
| `R/03_climate.R` | `s03_climate.py` | Pirson korrelyatsiyasi, tuman+yil ichki tahlil |
| `R/04_robustness.R` | `s04_robustness.py` | Jackknife, permutatsiya, qisman korrelyatsiya |
| `R/05_evt.R` | `s05_evt.py` | Gumbel, GEV, profil ehtimollik, bootstrap, VaR/CVaR |
| `R/06_product.R` | `s06_product.py` | To'lov, HE, bazis riski, tarif zonalari |
| `R/07_tables.R` | `s07_tables.py` | Panel FE, ta'sir kattaligi, 15 yillik oyna, EQN tariflari, tumanlar, issiqlik indeksi, shartli ozod |
| `run_all.R` | `run_all.py` | Butun zanjir; 285 ko'rsatkich `outputs/r_results.csv` ga yoziladi |
| `compare.R` | — | R va Python natijalarini kalitma-kalit solishtiradi |
| `reference/python_v1.1.0.csv` | — | Python v1.1.0 natijalari (285 ko'rsatkich) |
| `reference/make_reference.py` | — | Shu faylni Python'da qayta yaratish skripti |
| `data/` | `data/` | v1.1.0 dagi ma'lumotlarning aynan nusxasi |

## Qaysi vosita nimaga mos

| Hisob | Python | R |
|---|---|---|
| Trend, korrelyatsiya | `scipy.stats.linregress`, `pearsonr` | `lm`, `cor.test` |
| Panel FE, klaster SE | `statsmodels` (`cov_type="cluster"`, z-test) | `lm` + `sandwich::vcovCL(type = "HC1")` + `lmtest::coeftest(df = Inf)` |
| Gumbel va GEV (MLE) | `scipy.stats.gumbel_r.fit`, `genextreme.fit` | `evd::fgev` |
| Profil ehtimollik, LR test | o'z implementatsiyasi (`gev_profile_likelihood`) | `evd::profile` + `confint`, log-haqiqatga o'xshashliklar farqi |
| CVaR, EQN tariflari | `scipy.integrate.quad` | `integrate` |
| Protsentillar | `numpy.percentile` (linear) | `quantile(type = 7)` |
| Tasodifiy sonlar | NumPy PCG64, `seed = 42` | Mersenne-Twister, `set.seed(42)` |

## Solishtirish natijasi

`compare.R` har bir kalit uchun Python ko'rsatgan aniqlikni (`round(..., d)`) hisobga oladi.

| Holat | Soni | Ma'nosi |
|---|---|---|
| MOS | 273 | Deterministik natija Python bilan bir xil |
| MK: MOS | 9 | Bootstrap va permutatsiya — Monte-Karlo xatosi doirasida mos |
| OPTIMIZATOR | 3 | Erkin GEV bahosi — pastga qarang |
| FARQ | 0 | — |

Asosiy ko'rsatkichlar:

| Ko'rsatkich | Python | R |
|---|---|---|
| Panel kuzatuvlari | 254 | 254 |
| Panel FE (tuman+yil, 2018-yilsiz): β; p | −0,620; 0,013 | −0,6198; 0,0126 |
| Ichki R² | 0,0458 | 0,0458 |
| Chang gipotezasi: r; 2018-yilsiz r (p) | −0,820; −0,550 (0,100) | −0,8195; −0,5497 (0,0998) |
| Gumbel μ; σ (t) | −16 140; 28 024 | −16 140; 28 024 |
| Profil ehtimollik: ξ [95% CI] | 0,020 [−0,187; 0,373] | 0,0199 [−0,1865; 0,3727] |
| LR (Gumbel'ga nisbatan); p | 0,024; 0,876 | 0,0244; 0,8759 |
| 100 yillik daraja; CVaR(99%) | 19,0%; 23,7% | 18,96%; 23,68% |
| Bootstrap 95% CI (100 yil) | 10,0–28,5% | 10,2–28,5% |
| Bazis riski: zarar / qoplangan | 41 / 11 | 41 / 11 |
| HE: tuman / viloyat indeksi | 51,5% / 7,6% | 51,55% / 7,55% |
| Zonalar: barqaror / Boysun / yagona | 1,14 / 8,29 / 1,68% | 1,143 / 8,291 / 1,678% |

Maksimal haqiqatga o'xshashlik bahosiga asoslangan qiymatlar ikki mustaqil
optimizatorda 6 ta muhim raqamgacha mos keladi (nisbiy farq ~10⁻⁶).

## Tasodifiy sonlarga bog'liq natijalar

R va Python tasodifiy son generatorlari har xil: bir xil `seed = 42` bilan ham
ular boshqa ketma-ketlik beradi. Shuning uchun bootstrap oraliqlari va permutatsiya
p-qiymati raqamma-raqam emas, Monte-Karlo xatosi doirasida mos keladi:
oraliq chegaralari ±0,5 foiz punkti, permutatsiya p-qiymati ±0,003
(Python 0,0043, R 0,0047).

## Erkin GEV bahosi — muhim kuzatuv

| | ξ | σ | log-L |
|---|---|---|---|
| Python, `scipy.stats.genextreme.fit` | +1,141 | 0,2 t | −393,45 |
| R, `evd::fgev` | +0,020 | 27 895 t | −224,51 |
| Profil ehtimollik maksimumi (Python va R) | +0,020 | — | −224,51 |

`scipy` ning standart GEV moslashtirish funksiyasi bu qatorda haqiqiy maksimumni
topmaydi: u σ ≈ 0,2 t bo'lgan, taqsimotning quyi chegarasi eng kichik kuzatuvga
yopishgan nuqtada to'xtaydi. Bu nuqtaning log-haqiqatga o'xshashligi haqiqiy
maksimumdan 169 birlik past. Demak, ξ = +1,14 va bootstrap oralig'i
[−9,6; +8,6] parametr noaniqligini emas, optimizator nosozligini aks ettiradi.

R'da xuddi shu bootstrap (2000 takrorning hammasi muvaffaqiyatli yaqinlashgan) ξ uchun
[−0,98; +0,43] oraliqni beradi.

**Gumbel tanlovi bundan o'zgarmaydi** — u yanada mustahkam asosga ega bo'ladi:
GEV ning haqiqiy bahosi ξ̂ = 0,020, profil 95% oralig'i [−0,187; 0,373] nolni
o'z ichiga oladi va LR test Gumbel modelini rad etmaydi (p = 0,876). O'zgarishi
kerak bo'lgan narsa — asoslash: "ξ aniqlanmaydi, CI ≈ [−9,6; +8,6]" o'rniga
"ξ̂ ≈ 0 va Gumbel'dan statistik farq qilmaydi".

## Qisman korrelyatsiya haqida eslatma

Python (`s04_robustness.partial_correlation`) qoldiqlar korrelyatsiyasining
p-qiymatini n − 2 erkinlik darajasi bilan hisoblaydi. Bitta o'zgaruvchi nazorat
qilinganda to'g'ri qiymat n − 3. Shu sababli `run_all.R` ikkalasini ham chiqaradi:
r = −0,847, p = 0,0010 (n − 2) va p = 0,0020 (n − 3). Xulosa o'zgarmaydi.

## Python mos yozuvlarini qayta yaratish

SACRA repozitoriysining `src/` papkasi ko'rsatiladi:

```
python reference/make_reference.py /yo'l/sacra-index-insurance/src
```
