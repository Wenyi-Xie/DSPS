#' Simulating time-to-event data with time-dependent covariates
#'
#' \code{sim_data} is used to simulate longitudinal time-to-event data.
#' Users are allowed to specify the sample size, event rate and number of time-independent noise variables.
#' Currently, only one time-dependent covariate can be incorporated. The inclusion of multiple
#' time-dependent covariates will be permitted in the future versions.
#' Users can also generate long dataset at user-specified time-points.
#'The event time can be generated from two different distributions, either
#' from a proportional hazard model or an accelerated failure time model.
#'
#'
#' @param n the sample size
#' @param er the event rate (range 0-1)
#' @param dist a character variable that specify the distribution of event time.
#' "PH" for proportional hazard model (Default), "AFT" for accelerated failure time model.
#' @param noise the number of time-independent noise variables.
#' @param time_points an optional numeric vector containing the user-specified time points.
#'
#' @returns a list of three objects relating simulated time-to-event data with time-dependent covariates.
#'
#' \itemize{
#' \item{\code{T_failure}} is a subject-level dataframe that contains the underlying non-censored event time
#' and observation status (1: non-censored; 0: censored).
#'
#' \item{\code{Z}} is the longitudinal dataset of time-dependent covariates at unique event time
#' or user-specified time points (when the argument \code{time_point} is supplied).
#' Data points are preserved even for subjects who are no longer at risk.
#' The last four columns, in order, are \strong{time}, \strong{event status} (1:event; 0:non-event), \strong{unique identifier} for subject
#' and \strong{at-risk status} (1:at risk; 0: not at risk). The rest leading columns are time-dependent covariates.
#'
#' \item{\code{T_90}} is the 90th percentile of event time.
#' }
#'
#' @examples
#' set.seed(1236)
#'
#' dat <- sim_data(n=100, er=0.6)
#'
#' dat <- sim_data(n=200, er=0.3, dist="AFT", noise=10)
#'
#' dat <- sim_data(n=100, er=0.6, time_points = seq(from=0.1, to=1, by=0.1))
#'
#' @importFrom MASS mvrnorm
#' @import tidyverse
#' @import Matrix
#' @import dplyr
#'
#' @export
sim_data <- function(n, er, dist="PH", noise=0, time_points=NULL){

  n_all <- max(2000, n) # Total training sample size

  if(dist=="PH"){

    pho <- 0.5 # Pairwise correlation

    # Create covariates Covariance matrix and generate covariates Z
    Corr <- sapply(1:5, FUN = function(x) pho^{abs(x-1:5)})
    Sigma <- Corr*0.5^2
    Z <- mvrnorm(n=n_all, rep(0, 5), Sigma)

    # Covariate coefficients prespecified
    beta <- matrix(c(2, -1.6, 1.2, -0.8, 0.4)*1, ncol = 1)
    xgamma <- matrix(c(-0.4, 0.8, -1.2, 1.6, -2)*0, ncol=1)


    # Covariate coefficients for time-varying covariate
    betat <- 0.3

    sub_ee <- mvrnorm(n = n_all, mu = c(0, 1), Sigma = matrix(c(0.5^2, 0.5*0.5*0.2, 0.5*0.5*0.2, 0.5^2), nrow = 2))
    sub_effect <- sub_ee[, 1]
    sub_slope <- sub_ee[, 2]
    sub_slope <- ifelse(sub_slope<0, 0.1, sub_slope)
    # Generate true failure time by COX model
    U <- runif(n_all)

    # lambda <- 0.05 # Exponential parameter
    lambda <- 0.05
    T_failure <- log(1-(log(U)*sub_slope*betat)/(lambda*exp(Z %*% (beta + xgamma * betat)) * exp(sub_effect* betat)))/(sub_slope*betat)
  }else if(dist == "AFT"){

    pho <- 0.65 # Pairwise correlation

    beta_0 <- -2
    # Create covariates Covariance matrix and generate covariates Z
    Corr <- sapply(1:5, FUN = function(x) pho^{abs(x-1:5)})
    Sigma <- Corr*1^2
    Z <- mvrnorm(n=n_all, rep(0, 5), Sigma)

    # Covariate coefficients prespecified
    beta <- matrix(c(2, -1.6, 1.2, -0.8, 0.4)*0, ncol = 1)
    xgamma <- matrix(c(-0.4, 0.8, -1.2, 1.6, -2)*0, ncol=1)
    betaint <- matrix(c(0, -0.1, 0, 0.1, 0)*0, ncol = 1)

    # Covariate coefficients for time-varying covariate
    betat <- 0.8

    sub_ee <- mvrnorm(n = n_all, mu = c(0, 1), Sigma = matrix(c(0.5^2, 0.5*0.5*0.2, 0.5*0.5*0.2, 0.5^2), nrow = 2))
    sub_effect <- sub_ee[, 1]
    sub_slope <- sub_ee[, 2]
    sub_slope <- ifelse(sub_slope<0, 0.01, sub_slope)


    # Z_int <- cbind(Z[,1] * Z[, 1],
    #                Z[, 2] * Z[, 5],
    #                ifelse(Z[,4]^3 * sign(Z[,4]^3) <=2, Z[,4]^3, sign(Z[,4]^3)*2),
    #                abs(Z[, 1] + Z[, 2] + Z[, 3] + Z[, 4] + Z[, 5])
    #                )
    # Z_int_beta <- matrix(c(-0.25, 0.25, 0.5, -0.5)*1, ncol = 1)

    Z_int <- cbind(Z[,1] * Z[, 1],
                   Z[,2] * Z[, 5],
                   abs(Z[, 1] + Z[, 2] + Z[, 3] + Z[, 4] + Z[, 5]),
                   Z[,3] * Z[, 4])
    # Z_int_beta <- matrix(c(-0.5, -0.25, -0.8, -0.25)*1, ncol = 1)
    Z_int_beta <- matrix(c(-0.5, -0, -0.8, -0)*1, ncol = 1)

    #
    temp_slope <- sub_slope*(betat + Z %*% betaint)
    # Generate true failure time by COX model
    U <- rnorm(n_all)

    T_failure <- log(1+(exp(U)*temp_slope)/(exp(beta_0 + Z %*% beta + Z_int %*% Z_int_beta +  Z %*% betaint * sub_effect + sub_effect * betat)))/(temp_slope)

  }

  # hist(T_failure)
  # hist(T_failure1)
  # hist(T_failure2)
  #
  # #T_failure <- -log(U)/lambda*exp(-Z %*% beta)
  # quantile(T_failure, probs=c(0.9))
  # quantile(T_failure[grp == 1], probs=c(0.9))
  # quantile(T_failure[grp == 0], probs=c(0.9))


  T_90 <- quantile(T_failure, probs=c(0.9)) # 90th quantile --- censoring threshold
  T_90
  # hist(T_failure)

  # Generate censoring time, possibly two ways for doing this
  # Case one: AFT model
  # Case two: proportional hazard model

  Z <- Z[1:n, ]
  T_failure <- T_failure[1:n]
  # Censoring time Version 1 -- AFT model

  if(dist == "AFT"){
    # Z_int <- abs(Z[, 1] + Z[, 2] + Z[, 3] + Z[, 4] + Z[, 5])
    ########################################################################################.
    # Below are code that tests the correspondence of a and censoring rate ################.
    censor_ratio_1 <- function(x){
      T_cen_1 <- exp(sapply(Z %*% beta_c1 + x, FUN=function(x) rnorm(mean=x, sd=0.5, n=1)))
      # T_cen_1 <- exp(sapply(Z_int %*% beta_c1_int + x, FUN=function(x) rnorm(mean=x, sd=0.5, n=1)))
      ind_c1_ratio <- mean(I(T_failure < T_cen_1 & T_failure < T_90))
      return(ind_c1_ratio)
    }


    # Set event rate
    event_rate <- er
    beta_c1 <- matrix(c(1, 1, 1, -1, -1) * 0.5, ncol=1)
    # beta_c1_int <- matrix(0.2, ncol=1)
    c_ratio <- sapply(seq(0, 3, 0.01), censor_ratio_1)
    a <- seq(0, 3, 0.01)[rank(abs(c_ratio - event_rate), ties.method = "first") == 1]

    # T_cen <- exp(sapply(Z %*% beta_c1 + a, FUN=function(x) rnorm(mean=x, sd=0.5, n=1)))
    T_cen <- exp(sapply(Z %*% beta_c1 + a, FUN=function(x) rnorm(mean=x, sd=0.5, n=1)))

  }else if(dist == "PH"){

    ########################################################################################.
    # Below are code that tests the correspondence of a and censoring rate ################.
    censor_ratio_1 <- function(x){
      T_cen_1 <- exp(sapply(Z %*% beta_c1 + x, FUN=function(x) rnorm(mean=x, sd=0.5, n=1)))
      #T_cen_1 <- exp(sapply(Z_int %*% beta_c1_int + x, FUN=function(x) rnorm(mean=x, sd=0.5, n=1)))
      ind_c1_ratio <- mean(I(T_failure < T_cen_1 & T_failure < T_90))
      return(ind_c1_ratio)
    }


    # Set event rate
    event_rate <- er
    beta_c1 <- matrix(c(1, 1, 1, 1, 1) * 0.5, ncol=1)
    # beta_c1_int <- matrix(0.2, ncol=1)
    c_ratio <- sapply(seq(0, 3, 0.01), censor_ratio_1)
    a <- seq(0, 3, 0.01)[rank(abs(c_ratio - event_rate), ties.method = "first") == 1]

    T_cen <- exp(sapply(Z %*% beta_c1 + a, FUN=function(x) rnorm(mean=x, sd=0.5, n=1)))
    # T_cen <- exp(sapply(Z_int %*% beta_c1_int + a, FUN=function(x) rnorm(mean=x, sd=0.5, n=1)))
  }


  # Case 1
  # Observed time
  T_obs <- pmin(T_90, T_cen, T_failure) # For case two, switch T_cen_1 to T_cen_2
  # hist(T_obs)

  # Event indicator
  event_ind_obs <- 1*I(T_failure<T_cen & T_failure < T_90) # For case two, switch T_cen_1 to T_cen_2
  mean(event_ind_obs)



  if(noise > 0){
    Z_nuisance <- mvrnorm(n=n, rep(0, noise), diag(x=1, nrow=noise))
    Z <- cbind(Z, Z_nuisance)
  }


  uniq_t <- sort(unique(T_obs[event_ind_obs == 1]))

  kt <- sub_slope[1:n]
  sub_effect = sub_effect[1:n]


  if(!is.null(time_points)){
    uniq_t <- time_points
  }

  long_dat <- function(i){
    ## Keep only observed records, need to impute all time records when predicting
    # at_risk_n <- sum(T_obs[i] >=uniq_t)
    # if(at_risk_n > 0){
    #   Z_long <- t(replicate(at_risk_n, Z[i, ]))
    #   uniq_t_i <- uniq_t[T_obs[i] >=uniq_t]
    #
    #   Z_long <- as.data.frame(Z_long)
    #   Z_long <- cbind(Z_long,
    #                   uniq_t_i * kt[i] + sub_effect[i],
    #                   time = uniq_t_i,
    #                   event =ifelse(T_obs[i] == uniq_t_i, event_ind_obs[i], 0),
    #                   subjectid = paste0("S", i))
    #   return(Z_long)
    # }else{
    #   return(NULL)
    # }

    ## Keep all time records
    at_risk_n <- length(uniq_t)

    Z_long <- t(replicate(at_risk_n, Z[i, ]))

    Z_long <- as.data.frame(Z_long)
    Z_long <- cbind(Z_long,
                    uniq_t * kt[i] + sub_effect[i],
                    time = uniq_t,
                    event = ifelse(T_obs[i] == uniq_t, event_ind_obs[i], 0),
                    subjectid = paste0("S", i),
                    at_risk = 1*(T_obs[i] >=uniq_t))
    return(Z_long)

  }
  Z_long <- do.call("bind_rows", lapply(seq(n), long_dat))
  colnames(Z_long)[1:(length(colnames(Z_long)) - 4)] <- paste0("X", seq(length(colnames(Z_long)) - 4))




  base_data <- list(T_failure = data.frame(T_failure=ifelse(T_failure<T_90, T_failure, T_90),
                                           event_ind_obs=event_ind_obs,
                                           subjectid=paste0("S", seq(n))),
                    Z=Z_long,
                    T_90 = T_90)


  return(base_data)
}
