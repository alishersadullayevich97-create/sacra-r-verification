# 1-bosqich. Ma'lumotni yuklash va sifat qoidalarini qo'llash.

read_utf8 <- function(path) {
  read.csv(path, fileEncoding = "UTF-8", stringsAsFactors = FALSE, check.names = FALSE)
}

load_province <- function() {
  # Viloyat darajasidagi qator (2006-yil chiqarilgan).
  w <- read_utf8(WHEAT)
  v <- w[w$district == PROVINCE_ROW & w$year != EXCLUDE_YEAR, ]
  v <- v[order(v$year), ]
  rownames(v) <- NULL
  v
}

load_districts <- function() {
  # Tumanlar paneli. Qoidalar (dissertatsiya 2.1.3):
  #   - 2006-yil butunlay chiqariladi;
  #   - Bandixon 2011-2019 yillarda Qiziriq tarkibida bo'lgani uchun chiqariladi;
  #   - 2010-yil Bandixon va Qiziriq kuzatuvlari chiqariladi (v1.1.0);
  #   - Qiziriq 2024-yil kuzatuvi anomaliya sifatida chiqariladi;
  #   - viloyat va shahar qatorlari tahlilga kirmaydi.
  w <- read_utf8(WHEAT)
  w <- w[w$year != EXCLUDE_YEAR, ]
  w <- w[!(w$district == "Bandixon" &
           w$year >= BANDIXON_MERGED[1] & w$year <= BANDIXON_MERGED[2]), ]
  w <- w[!w$district %in% c(PROVINCE_ROW, CITY_ROW), ]
  for (i in seq_len(nrow(DROP_OBS))) {
    hit <- w$district == DROP_OBS$district[i] & w$year == DROP_OBS$year[i]
    w$yield_ts_ga[hit] <- NA
  }
  w <- w[!is.na(w$yield_ts_ga), ]
  rownames(w) <- NULL
  w
}

load_climate <- function() read_utf8(CLIMATE)
load_station <- function() read_utf8(STATION)
