estimate_simes_power_vec <- function(h, s, R, nreps_h, m, alpha) {
  theta <- numeric(m)
  if (s > 0) theta[1:s] <- h / s
  
  Z <- matrix(rnorm(m * nreps_h), nrow = m, ncol = nreps_h)
  Xmat <- matrix(theta, nrow = m, ncol = nreps_h) + t(R) %*% Z
  P <- pnorm(Xmat, lower.tail = FALSE)
  
  simes_rej <- apply(P, 2, function(colp) {
    sp <- sort(colp)
    min(sp * (m / (1:m))) < alpha
  })
  mean(simes_rej)
}

find_h_for_target_power_vec <- function(s, R, m, target_power = 0.95, nreps_h = 500,
                                        tol_power = 0.01, max_iter = 30, alpha = 0.05) {
  h_low <- 0
  power_low <- estimate_simes_power_vec(h_low, s, R, nreps_h, m, alpha)
  if (power_low >= target_power) {
    return(0)
  }
  
  h_high <- 0.1
  power_high <- estimate_simes_power_vec(h_high, s, R, nreps_h, m, alpha)
  iter <- 0
  
  while (power_high < target_power && iter < max_iter) {
    h_high <- h_high * 2
    if (h_high == 0) h_high <- 0.1
    power_high <- estimate_simes_power_vec(h_high, s, R, nreps_h, m, alpha)
    iter <- iter + 1
  }
  
  if (power_high < target_power) {
    warning(sprintf("Couldn't reach target power=%.3f for s=%d; returning h_high=%.4g (power=%.3f)",
                    target_power, s, h_high, power_high))
    return(h_high)
  }
  
  for (it in seq_len(50)) {
    h_mid <- (h_low + h_high) / 2
    power_mid <- estimate_simes_power_vec(h_mid, s, R, nreps_h, m, alpha)
    if (abs(power_mid - target_power) <= tol_power) {
      return(h_mid)
    }
    if (power_mid < target_power) {
      h_low <- h_mid
      power_low <- power_mid
    } else {
      h_high <- h_mid
      power_high <- power_mid
    }
  }
  
  (h_low + h_high) / 2
}