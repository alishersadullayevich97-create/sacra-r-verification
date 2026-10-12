# R natijalarini (outputs/r_results.csv) Python SACRA v1.1.0 natijalari
# (reference/python_v1.1.0.csv) bilan solishtiradi.
#
# Ishga tushirish: avval run_all.R, keyin  Rscript compare.R  yoki  source("compare.R")
#
# Holatlar:
#   MOS         - deterministik qiymat, Python ko'rsatgan aniqlik doirasida bir xil
#   FARQ        - deterministik qiymat farq qiladi (tekshirish kerak)
#   MK: MOS     - tasodifiy son generatoriga bog'liq (bootstrap, permutatsiya);
#                 farq Monte-Karlo xatosi doirasida
#   MK: FARQ    - Monte-Karlo chegarasidan tashqarida
#   OPTIMIZATOR - natija optimallashtirish algoritmiga bog'liq (erkin GEV); faqat ma'lumot uchun

find_root_cmp <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f) == 1) return(dirname(normalizePath(f)))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(dirname(normalizePath(of)))
  }
  normalizePath(".")
}
ROOT_CMP <- find_root_cmp()

py <- read.csv(file.path(ROOT_CMP, "reference", "python_v1.1.0.csv"), colClasses = "character", fileEncoding = "UTF-8")
rr <- read.csv(file.path(ROOT_CMP, "outputs", "r_results.csv"), colClasses = "character", fileEncoding = "UTF-8")
names(py)[2] <- "python"; names(rr)[2] <- "r"
cmp <- merge(py, rr, by = "key", all = TRUE)

# Python kodidagi round(..., d) bo'yicha ko'rsatilgan aniqlik (o'nli xonalar soni).
decimals <- list(
  "^pearson\\..*\\.r$" = 3, "^pearson\\..*\\.p$" = 4,
  "^dust\\.jackknife\\..*\\.r$" = 3, "^dust\\.jackknife_worst_r$" = 3, "^dust\\.jackknife_worst_p$" = 4,
  "^gev\\.profile\\.(xi|ci_low|ci_high|LR|p)$" = 3, "^gev\\.profile\\.ll_" = 2,
  "^rl\\..*\\.t$" = 0, "^rl\\..*\\.pct$" = 1, "^c?var[0-9]+\\.pct$" = 1, "^window2010\\.(rl|cvar)" = 1,
  "^index\\..*\\.sof_tarif_pct$" = 2, "^index\\..*\\.tolov_chastotasi_pct$" = 0,
  "^index\\..*\\.HE_pct$" = 1, "^index\\..*\\.korrelyatsiya$" = 2,
  "^zone\\..*\\.(sof|yuklamali|subsidiya|fermer)$" = 2, "^zone\\..*\\.HE_pct$" = 1,
  "^zone\\..*\\.chastota_pct$" = 0, "^zone\\..*\\.tumanlar$" = 0,
  "^fe\\." = 3,
  "^effect\\.(within_sd_kun|tasir_pct)$" = 2, "^effect\\.(r|R2)$" = 3,
  "^event2018\\.kvantil_pct$" = 2, "^event2018\\.davr_yil$" = 0, "^event2018\\.P_n_yilda_pct$" = 1,
  "^evt_tariff\\." = 2, "^province\\.yield_anom_sd$" = 1,
  "^district\\..*\\.(sd_pct|HE_tuman|HE_viloyat)$" = 1, "^district\\..*\\.r_viloyat$" = 2,
  "^district\\..*\\.sof_tarif$" = 2, "^district\\..*\\.(chastota_pct|n)$" = 0,
  "^heat_index\\..*\\.sof_tarif$" = 2, "^heat_index\\..*\\.chastota_pct$" = 0,
  "^heat_index\\..*\\.HE_pct$" = 1, "^heat_index\\..*\\.korrelyatsiya$" = 2,
  "^deductible\\..*\\.sof_tarif_pct$" = 2, "^deductible\\..*\\.tolov_chastotasi_pct$" = 0,
  "^deductible\\..*\\.HE_pct$" = 1, "^deductible\\..*\\.korrelyatsiya$" = 2,
  "\\.n$|\\.n_|_events$|\\.covered$|\\.uncovered$|\\.false_payouts$|worst_year$" = 0
)
mc_tol <- list("^boot\\." = 0.5,                  # foiz punkti; bootstrap percentil oralig'i
               "^dust\\.permutation_p$" = 0.003)   # ~3 Monte-Karlo standart xatosi
optimizer_keys <- "^gev\\.free_"

lookup <- function(key, table) {
  for (pat in names(table)) if (grepl(pat, key)) return(table[[pat]])
  NULL
}

judge <- function(key, p, r) {
  if (is.na(p) || is.na(r)) return(c("YO'Q", NA))
  pn <- suppressWarnings(as.numeric(p)); rn <- suppressWarnings(as.numeric(r))
  if (is.na(pn) || is.na(rn)) return(c(if (identical(p, r)) "MOS" else "FARQ", NA))
  d <- rn - pn
  if (grepl(optimizer_keys, key)) return(c("OPTIMIZATOR", d))
  mc <- lookup(key, mc_tol)
  if (!is.null(mc)) return(c(if (abs(d) <= mc) "MK: MOS" else "MK: FARQ", d))
  dec <- lookup(key, decimals)
  tol <- if (!is.null(dec)) 0.6 * 10^(-dec) else max(1e-4 * abs(pn), 1e-9)   # yaxlitlanmagan: 4 muhim xona
  c(if (abs(d) <= tol) "MOS" else "FARQ", d)
}

j <- t(mapply(judge, cmp$key, cmp$python, cmp$r))
cmp$farq  <- suppressWarnings(as.numeric(j[, 2]))
cmp$holat <- j[, 1]
write.csv(cmp, file.path(ROOT_CMP, "outputs", "comparison.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("\nR va Python (SACRA v1.1.0) natijalarini solishtirish\n", strrep("-", 60), "\n", sep = "")
print(table(holat = cmp$holat))
other <- cmp[cmp$holat != "MOS", c("key", "python", "r", "holat")]
if (nrow(other)) {
  cat("\nMOS'dan boshqa holatdagi ko'rsatkichlar:\n")
  other$python <- suppressWarnings(formatC(as.numeric(other$python), digits = 6, format = "g"))
  other$r      <- suppressWarnings(formatC(as.numeric(other$r),      digits = 6, format = "g"))
  print(other, row.names = FALSE)
}

headline <- c("panel.n_obs", "fe.twoway_fe.beta_2018siz", "fe.twoway_fe.p_2018siz", "within.R2",
              "dust.r", "dust.jackknife_worst_r", "dust.jackknife_worst_p",
              "gumbel.mu", "gumbel.sigma", "gev.profile.xi", "gev.profile.ci_low", "gev.profile.ci_high",
              "gev.profile.LR", "gev.profile.p", "rl.T100.pct", "cvar99.pct",
              "basis.loss_events", "basis.covered", "index.district.HE_pct", "index.province.HE_pct",
              "zone.stable.sof", "zone.high.sof", "zone.single.sof", "zone.high.fermer", "zone.high_risk")
h <- cmp[match(headline, cmp$key), c("key", "python", "r", "holat")]
fmt <- function(x) {
  n <- suppressWarnings(as.numeric(x))
  ifelse(is.na(n), x, suppressWarnings(formatC(n, digits = 6, format = "g")))
}
h$python <- fmt(h$python); h$r <- fmt(h$r)
cat("\nAsosiy ko'rsatkichlar:\n"); print(h, row.names = FALSE)

n_bad <- sum(cmp$holat %in% c("FARQ", "MK: FARQ", "YO'Q"))
cat(sprintf("\nJami %d ta ko'rsatkich. %s\n", nrow(cmp),
            if (n_bad == 0) "Deterministik natijalarning barchasi Python bilan mos keldi."
            else sprintf("%d ta ko'rsatkichni tekshirish kerak (outputs/comparison.csv).", n_bad)))
