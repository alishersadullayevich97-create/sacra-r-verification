# 4-bosqich. Iqlim omillari va hosildorlik anomaliyasi o'rtasidagi bog'liqlik.

pearson_table <- function(df, y_col, x_cols) {
  do.call(rbind, lapply(x_cols, function(x) {
    sub <- df[complete.cases(df[, c(y_col, x)]), ]
    ct  <- cor.test(sub[[y_col]], sub[[x]], method = "pearson")
    data.frame(omil = x, r = unname(ct$estimate), p = ct$p.value, n = nrow(sub),
               ahamiyatli = ct$p.value < 0.05, stringsAsFactors = FALSE)
  }))
}

two_way_within <- function(df, y_col, x_col, unit = "district", time = "year") {
  # Tuman o'rtachasi, so'ng yil o'rtachasi ayiriladi (bir martalik ikki bosqichli
  # markazlashtirish - Python s03_climate.two_way_within bilan aynan bir xil).
  d  <- df[complete.cases(df[, c(y_col, x_col)]), ]
  dm <- function(v) {
    u <- v - ave(v, d[[unit]])
    u - ave(u, d[[time]])
  }
  ct <- cor.test(dm(d[[y_col]]), dm(d[[x_col]]))
  r  <- unname(ct$estimate)
  list(r = r, p = ct$p.value, r_squared = r^2, n = nrow(d))
}
