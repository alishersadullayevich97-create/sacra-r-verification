# 5-8-bosqichlar. Ekstremal qiymatlar nazariyasi, bootstrap va risk o'lchovlari
# (dissertatsiya 2.3-bo'limi).
#
# Baholash evd paketi bilan bajariladi (fgev - maksimal haqiqatga o'xshashlik).
# Yo'qotishlar tonnada (o'n minglab) bo'lgani uchun optimallashtirish ichki
# ravishda x / s masshtabda olib boriladi (s - qatorning SD si). Lokatsiya-masshtab
# MLE bu o'zgarishga ekvivariant: mu va sigma s ga ko'paytirib qaytariladi,
# log-haqiqatga o'xshashlikdan n*log(s) ayiriladi. Natija asl birlikdagi MLE bilan
# aynan bir xil, faqat optimizator barqarorroq ishlaydi.

loss_series <- function(anomaly_abs) -anomaly_abs   # musbat qiymat = yo'qotish

# optim (BFGS) uchun qat'iy yaqinlashish mezoni: standart reltol = 1e-8 yetarli, ammo
# Python (scipy) aniq MLE bilan tonnagacha mos kelishi uchun qattiqroq qo'yiladi.
OPTIM_CONTROL <- list(reltol = 1e-12, maxit = 10000)

fit_gumbel <- function(x, s = sd(x)) {
  f <- suppressWarnings(evd::fgev(x / s, shape = 0, std.err = FALSE, control = OPTIM_CONTROL))
  if (f$convergence != "successful") stop("Gumbel MLE yaqinlashmadi")
  list(mu = f$estimate[["loc"]] * s, sigma = f$estimate[["scale"]] * s,
       loglik = as.numeric(logLik(f)) - length(x) * log(s), fit = f, s = s)
}

fit_gev <- function(x, s = sd(x), std.err = FALSE) {
  f <- suppressWarnings(evd::fgev(x / s, std.err = std.err, control = OPTIM_CONTROL))
  if (f$convergence != "successful") stop("GEV MLE yaqinlashmadi")
  list(xi = f$estimate[["shape"]], mu = f$estimate[["loc"]] * s,
       sigma = f$estimate[["scale"]] * s,
       loglik = as.numeric(logLik(f)) - length(x) * log(s), fit = f, s = s)
}

check_gev_stability <- function(losses, n_boot = N_BOOTSTRAP, seed = RANDOM_SEED) {
  # Erkin GEV bahosi va shakl parametrining bootstrap 95% oralig'i.
  s    <- sd(losses)
  base <- fit_gev(losses, s)
  set.seed(seed)
  xis <- replicate(n_boot, {
    b <- sample(losses, length(losses), replace = TRUE)
    tryCatch(fit_gev(b, s)$xi, error = function(e) NA_real_)
  })
  ok <- xis[is.finite(xis)]
  ci <- quantile(ok, c(0.025, 0.975), type = 7, names = FALSE)
  list(xi = base$xi, sigma = base$sigma, ci_low = ci[1], ci_high = ci[2],
       identified = !(ci[1] < 0 && 0 < ci[2]), n_ok = length(ok), n_boot = n_boot)
}

gev_profile_likelihood <- function(losses, level = 0.95) {
  # Shakl parametrining profil ehtimolligi (evd::profile + confint) va
  # Gumbel (xi = 0) modeliga nisbatan ehtimollik nisbati testi (evd::anova).
  s   <- sd(losses)
  gev <- fit_gev(losses, s, std.err = TRUE)   # profile() standart xatoni talab qiladi
  gum <- fit_gumbel(losses, s)
  invisible(capture.output(                    # evd o'z xabarini print() bilan chiqaradi
    pr <- profile(gev$fit, which = "shape", conf = 0.999)))
  invisible(capture.output(ci <- confint(pr, parm = "shape", level = level)))
  lr  <- max(0, 2 * (gev$loglik - gum$loglik))
  list(xi = gev$xi, ci_low = ci[1, 1], ci_high = ci[1, 2],
       ll_gev = gev$loglik, ll_gumbel = gum$loglik,
       LR = lr, p = pchisq(lr, df = 1, lower.tail = FALSE))
}

return_levels <- function(mu, sigma, periods = c(10, 20, 50, 100), mean_production = NULL) {
  # T yillik qaytarilish darajasi: x_T = mu - sigma * ln(-ln(1 - 1/T))
  q <- evd::qgumbel(1 - 1 / periods, loc = mu, scale = sigma)
  out <- data.frame(qaytarilish_davri = periods, yoqotish_t = q)
  if (!is.null(mean_production)) out$yoqotish_pct <- q / mean_production * 100
  out
}

bootstrap_return_levels <- function(losses, periods = c(10, 20, 50, 100), mean_production = NULL,
                                    n_boot = N_BOOTSTRAP, seed = RANDOM_SEED) {
  # Qaytarilish darajalarining 95% bootstrap oraliqlari (percentil usuli).
  s     <- sd(losses)
  scale <- if (is.null(mean_production)) 1 else 100 / mean_production
  set.seed(seed)
  acc <- t(replicate(n_boot, {
    b <- sample(losses, length(losses), replace = TRUE)
    g <- tryCatch(fit_gumbel(b, s), error = function(e) NULL)
    if (is.null(g)) rep(NA_real_, length(periods))
    else evd::qgumbel(1 - 1 / periods, g$mu, g$sigma) * scale
  }))
  base <- fit_gumbel(losses, s)
  pt   <- evd::qgumbel(1 - 1 / periods, base$mu, base$sigma) * scale
  ci   <- apply(acc, 2, function(a) quantile(a[is.finite(a)], c(0.025, 0.975), type = 7, names = FALSE))
  data.frame(qaytarilish_davri = periods, nuqtaviy = pt, quyi_2.5 = ci[1, ], yuqori_97.5 = ci[2, ],
             check.names = FALSE)
}

var_cvar <- function(mu, sigma, alphas = c(0.95, 0.99), mean_production = NULL) {
  # VaR - kvantil; CVaR = 1/(1-a) * int_a^0.999999 F^-1(q) dq
  # (yuqori chegara Python bilan bir xil: 0.999999).
  rows <- lapply(alphas, function(a) {
    var  <- evd::qgumbel(a, mu, sigma)
    cv   <- integrate(function(q) evd::qgumbel(q, mu, sigma), a, 0.999999,
                      subdivisions = 200L, rel.tol = 1e-10)$value
    cvar <- cv / (1 - a)
    row  <- data.frame(ishonch = a, VaR_t = var, CVaR_t = cvar)
    if (!is.null(mean_production)) {
      row$VaR_pct  <- var  / mean_production * 100
      row$CVaR_pct <- cvar / mean_production * 100
    }
    row
  })
  do.call(rbind, rows)
}
