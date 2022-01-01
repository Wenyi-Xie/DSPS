
mosek_solve <- function(x, y, Cn, kparam, weight, kernel){
  n_at_risk_i_j <- length(y)
  d <- matrix(data=1, nrow = n_at_risk_i_j)
  Beq <- matrix(rep(0, 1), ncol = 1)
  b <- weight * Cn

  prob <- list(sense="min")
  prob$c <- -d
  prob$A <- as(matrix(0, nrow=1, ncol=n_at_risk_i_j), "sparseMatrix")
  prob$bc <- t(cbind(blc=Beq, buc=Beq))
  prob$bx <- t(cbind(blx=rep(0, n_at_risk_i_j), bux = b))

  if(kernel == "semi_gaussian"){
    xinner <- xinner_kernel(x = x, y =x, kernel=kernel, kparam = kparam)
  }else if(kernel == "semi_linear"){
    xinner <- xinner_linear(x = x, y = x)
  }

  C_tv <- y * (xinner + 1)
  C_tv <- as(t(t(C_tv) * y), "dgTMatrix")


  prob$qobj$i <- C_tv@i[C_tv@i>=C_tv@j] + 1
  prob$qobj$j <- C_tv@j[C_tv@i>=C_tv@j] + 1
  prob$qobj$v <- C_tv@x[C_tv@i>=C_tv@j]  + 10^(-6)*I(C_tv@i == C_tv@j)[C_tv@i>=C_tv@j]
  prob$iparam <- list(LOG = 0)
  prob$dparam <- list(INTPNT_QO_TOL_DFEAS = 10^-5,
                      INTPNT_QO_TOL_PFEAS = 10^-5)

  cvhm <- mosek(prob)


  betahat <- cvhm$sol$itr$xx  * y
  alpha_hat <- sum(cvhm$sol$itr$xx * y)

  fit.mosek <- list(beta = matrix(betahat, ncol = 1), beta0 = alpha_hat)
  return(fit.mosek)
}
