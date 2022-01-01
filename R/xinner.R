Kmat <- function(x, 
                 y,  
                 kernel = "gaussian",  
                 kparam = 1.0) {
  ncol = dim(y)[2]
  x <- matrix(x, ncol= ncol)
  
  if( kernel == "polynomial" ) {
    obj <- (x %*% t(y) + 1.0)^kparam
  } else if( kernel == "gaussian" ) {
    normx <- drop((x^2) %*% rep(1.0, ncol(x)))
    normy <- drop((y^2) %*% rep(1.0, ncol(y)))
    temp <- x %*% t(y)
    temp <- (-2.0 * temp + normx) + outer(rep(1.0, nrow(x)), normy, "*")
    obj <- exp(-temp / (2*kparam^2))
  } else if(kernel == "semi_gaussian"){
    
    ytemp = y[, -c(ncol)]
    ytime = y[, ncol] # add time
    
    xtemp = matrix(x[, 1:(ncol-1)], ncol=ncol-1)
    xtime = x[, ncol] # add time
    
    time_out <- outer(rep(1.0, length(xtime)), ytime)
    time_kernel <- ifelse(time_out == xtime, 1, 0)
    
    normx <- drop((xtemp^2) %*% rep(1.0, ncol(xtemp)))
    normy <- drop((ytemp^2) %*% rep(1.0, ncol(ytemp)))
    temp <- xtemp %*% t(ytemp)
    
    temp <- (-2.0 * temp + normx) + outer(rep(1.0, nrow(x)), normy, "*")
    obj <- exp(-temp / (2*kparam^2)) + time_kernel
  }
  
  return(obj)
  
}

xinner_kernel <- function(x, y, kernel, kparam){

  K.train <- Kmat(x = x, 
                  y = y, 
                  kernel = kernel, 
                  kparam = kparam)
  
  xinner <- K.train
  
  return(xinner)
  
}

xinner_linear <- function(x, y){
  ncol = dim(y)[2]
  ytemp = y[, -c(ncol)]
  ytime = y[, ncol] # add time
  
  x <- matrix(x, ncol= ncol)
  xtemp = matrix(x[, 1:(ncol-1)], ncol=ncol-1)
  xtime = x[, ncol] # add time
  
  time_out <- outer(rep(1.0, length(xtime)), ytime)
  time_kernel <- ifelse(time_out == xtime, 1, 0)
  
  xinner <- xtemp %*% t(ytemp) + time_kernel

  return(xinner)
}
