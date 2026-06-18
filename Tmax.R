Tmax <- function(Y1, Y2, K = 4) {
  # Ordered categories assumed to be 1,2,3,4 (K=4 by default).
  # Works even if some categories are not observed.
  
  # drop NAs
  Y1 <- Y1[!is.na(Y1)]
  Y2 <- Y2[!is.na(Y2)]
  
  # basic checks
  if (length(Y1) == 0 || length(Y2) == 0) return(NA_real_)
  if (any(!Y1 %in% 1:K) || any(!Y2 %in% 1:K)) {
    stop("Y1 and Y2 must take values in {1,2,...,K}.")
  }
  
  # tabulate counts for all categories 1,...,K
  r1 <- tabulate(Y1, nbins = K)
  r2 <- tabulate(Y2, nbins = K)
  t  <- r1 + r2
  
  # keep only pooled-observed categories
  idx <- t > 0
  c   <- sum(idx)
  
  n1 <- sum(r1)
  n2 <- sum(r2)
  if (n1 == 0 || n2 == 0) return(NA_real_)
  
  # group indicator
  W <- rep(c(0, 1), times = c(n1, n2))
  
  # if only one category is observed overall
  if (c <= 1) return(0)
  
  # stochastic ordering check
  lhs <- cumsum(r2[idx]) / n2
  rhs <- cumsum(r1[idx]) / n1
  
  if (sum(lhs >= rhs) == c) {
    
    M <- matrix(0, nrow = c - 1, ncol = c)
    E <- matrix(0, nrow = c - 1, ncol = n1 + n2)
    F <- numeric(c - 1)
    
    for (i in 1:(c - 1)) {
      M[i, ] <- rep(c(0, 1), times = c(i, c - i))
      E[i, ] <- rep(c(M[i, ], M[i, ]),
                    times = c(r1[idx], r2[idx]))
      F[i]   <- suppressWarnings(cor(E[i, ], W, method = "pearson"))
    }
    
    T <- max(F, na.rm = TRUE)
    
  } else {
    
    # isotonic regression branch
    if (!exists("pava", mode = "function")) {
      stop("pava() not found. Please source pava.R or load a package providing pava().")
    }
    
    y   <- r2[idx] / t[idx]
    w   <- t[idx]
    iso <- pava(y, w)
    
    Yt  <- rep(c(iso, iso),
               times = c(r1[idx], r2[idx]))
    T   <- suppressWarnings(cor(Yt, W, method = "pearson"))
  }
  
  as.numeric(T)
}