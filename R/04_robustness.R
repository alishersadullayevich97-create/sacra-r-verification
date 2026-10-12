# Barqarorlik sinovlari (dissertatsiya 2.2.7 va B ilova): chang bo'roni gipotezasi.

jackknife_correlation <- function(x, y) {
  # Har bir kuzatuvni navbat bilan chiqarib, korrelyatsiyani qayta hisoblaydi.
  ct0 <- cor.test(x, y)
  r0  <- unname(ct0$estimate)
  jk  <- do.call(rbind, lapply(seq_along(x), function(i) {
    ct <- cor.test(x[-i], y[-i])
    data.frame(chiqarilgan_indeks = i, r = unname(ct$estimate), p = ct$p.value,
               ahamiyatli = ct$p.value < 0.05, r_ozgarishi = unname(ct$estimate) - r0)
  }))
  list(r0 = r0, p0 = ct0$p.value, table = jk)
}

permutation_test <- function(x, y, n_perm = 10000, seed = RANDOM_SEED) {
  # Taqsimot farazisiz p-qiymat. R generatori Python'nikidan farq qiladi,
  # shuning uchun natija Monte-Karlo xatosi doirasida mos keladi, raqamma-raqam emas.
  set.seed(seed)
  obs <- abs(cor(x, y))
  hits <- replicate(n_perm, abs(cor(x, sample(y))) >= obs)
  mean(hits)
}

partial_correlation <- function(x, y, z) {
  # z nazorat qilingan qisman korrelyatsiya (qoldiqlar korrelyatsiyasi).
  # p_python: Python bilan bir xil (erkinlik darajasi n-2);
  # p_dof:    nazorat o'zgaruvchisi hisobga olingan to'g'ri variant (n-3).
  rx <- resid(lm(x ~ z)); ry <- resid(lm(y ~ z))
  r  <- cor(rx, ry)
  n  <- length(x)
  t2 <- function(df) 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df)
  list(r = r, p_python = t2(n - 2), p_dof = t2(n - 3))
}
