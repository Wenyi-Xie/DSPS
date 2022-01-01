#' Fitting (Training) dynamic survival prediction
#'
#' \code{dsp_fit} is used to fit the dynamic survival prediction.
#' It allows for both linear and Gaussian kernel as well as two time-discrete kernel, semi-linear and semi-Gaussian kernel as described in the manuscript.
#' A list of estimated parameters are returned.
#'
#' @param data a data frame with time-dependent variables at unique event time points.
#' \code{time}, \code{event}, \code{subjectid} should be included after the leading columns of time-dependent variables (Order doesn't matter).
#' \code{time} is the ordered unique event time points. \code{event} is the survival status at \code{time} (1:event; 0:non-event). \code{subjectid} is the unique identifier for subject.
#'
#' @param kernel the kernel to be used; Four types of kernels are currently supported: "linear", "gaussian", "semi_linear", "semi_gaussian"
#'
#' @param sigma a positive number, the bandwidth parameter for "gaussian" and "semi_gaussian" kernel.
#' @param C a positive number, the penalty parameter that controls the penalty of mis-classification.
#'
#' @return A list of fitting results that can be passed to function \code{dsp_predict}.
#'
#' @examples
#'
#' data_list <- sim_data(n=100, er=0.6)
#'
#' # Keep only at-risk observations
#' train_data <- data_list$Z %>%
#'                 filter(at_risk == 1) %>%
#'                 select(-at_risk)
#'
#' fit <- dsp_fit(train_data,
#'                kernel="semi_linear",
#'                C=2^0)
#'
#' fit <- dsp_fit(train_data,
#'                kernel="semi_gaussian",
#'                sigma=2^0,
#'                C=2^0)
#'
#' @import Rmosek
#'
#' @export

dsp_fit <- function(data,
                    kernel,
                    C,
                    sigma=NULL){


  transformed_dat <- data_transform(data)

  # Standardize the long covariate matrix
  train_mean <- colMeans(transformed_dat$x)
  train_sd <- apply(transformed_dat$x, MARGIN = 2, sd)
  trainx <- sweep(sweep(transformed_dat$x, 2, train_mean), 2, train_sd, FUN="/")
  transformed_dat$x <- as.matrix(trainx)

  if(kernel=="semi_gaussian" |
     kernel=="semi_linear"){
    fit.dsp <- mosek_solve(x = transformed_dat$x,
                          y = transformed_dat$y,
                          Cn = C,
                          kparam = sigma,
                          weight = transformed_dat$weight,
                          kernel = kernel)
  }else if(kernel=="gaussian"){
    fit.dsp <- wsvm(x=transformed_dat$x,
                       y=as.factor(transformed_dat$y),
                       gamma=1/(2*sigma^2),
                       cost = C,
                       kernel = "radial",
                       weight = transformed_dat$weight,
                       tolerance = 0.0001)
  }else if(kernel=="linear"){
    fit.dsp <- wsvm(x=transformed_dat$x,
                       y=as.factor(transformed_dat$y),
                       # gamma=1/(2*kparam^2),
                       cost = C,
                       kernel = "linear",
                       weight = transformed_dat$weight)
  }

  fit.dsp$sigma = sigma
  fit.dsp$C = C
  fit.dsp$kernel = kernel
  fit.dsp$train_dat = transformed_dat
  fit.dsp$train_mean = train_mean
  fit.dsp$train_sd = train_sd
  return(fit.dsp)
}
