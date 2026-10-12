# Markaziy sozlamalar. Qiymatlar SACRA v1.1.0 dagi src/config.py bilan bir xil.
# Barcha modullar parametrlarni shu yerdan oladi.

DATA <- file.path(ROOT, "data")
OUT  <- file.path(ROOT, "outputs")
dir.create(OUT, showWarnings = FALSE)

# --- Ma'lumot fayllari ---
WHEAT   <- file.path(DATA, "wheat_panel.csv")
CLIMATE <- file.path(DATA, "climate_annual_era5.csv")
STATION <- file.path(DATA, "station_annual.csv")

# --- Ma'lumot sifati bo'yicha qarorlar (dissertatsiya 2.1.3-bo'limi) ---
EXCLUDE_YEAR    <- 2006                 # hisobot metodologiyasidagi uzilish
BANDIXON_MERGED <- c(2011, 2019)        # tuman Qiziriq tarkibida bo'lgan davr
PROVINCE_ROW    <- "Surxondaryo (viloyat)"
CITY_ROW        <- "Termiz shahri"
# Alohida chiqariladigan kuzatuvlar: Qiziriq 2024 anomaliyasi va 2010-yildagi
# Bandixon-Qiziriq o'tish kuzatuvlari (v1.1.0).
DROP_OBS <- data.frame(district = c("Qiziriq", "Bandixon", "Qiziriq"),
                       year     = c(2024,      2010,       2010),
                       stringsAsFactors = FALSE)

# --- Model parametrlari ---
MIN_YEARS_TREND <- 8      # trend baholash uchun minimal kuzatuv soni
N_BOOTSTRAP     <- 2000   # bootstrap takrorlari
RANDOM_SEED     <- 42     # R va Python tasodifiy son generatorlari har xil (README ga qarang)

# --- Sug'urta mahsuloti parametrlari (dissertatsiya 3.1-3.2 bo'limlari) ---
DEDUCTIBLE     <- 10.0    # shartli ozod, %
RISK_LOADING   <- 0.25    # risk yuklamasi
SUBSIDY_RATE   <- 0.50    # davlat subsidiyasi ulushi (VM 206-son qarori, 15-ilova)
SUBSIDY_CAP    <- NA      # NA - chegarasiz; avvalgi tartib (VM 373-son qarori) uchun 4.0
LOSS_THRESHOLD <- 10.0    # bazis riski tahlilida zarar chegarasi, %

# --- Tarif zonalari (dissertatsiya 3.2.3-bo'limi) ---
ZONE_SD_MIN    <- 14.0    # anomaliya SD > 14%
ZONE_R_MAX     <- 0.5     # va viloyat bilan korrelyatsiya < 0,5
HIGH_RISK_ZONE <- NULL    # NULL - mezon bo'yicha; tumanlar vektori berilsa, shu ishlatiladi
