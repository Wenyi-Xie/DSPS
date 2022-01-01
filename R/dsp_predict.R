#' Predicting the event time from DSP
#'
#' Predicted event time based on DSP method.
#'
#' @param fit is the return object from function \code{dsp_fit}
#' @param newdata a dataframe in which to look for time-dependent covariates with which to predict.
#' The last two columns, in order, are time (\code{time}) and unique identifier for subject (\code{subjectid}).
#' The order of the rest time-dependent covariates should match the order as in the training data used in function \code{dsp_fit}.
#' For each subject, time-dependent covariates should be supplied at each unique event time in the training data.
#' Missing values should be imputed before making prediction.
#'
#' @return a dataframe with predicted event time
#'
#' @examples
#'
#' ## Example 1 - With training data splitting
#'
#' set.seed(1236)
#' dat <- sim_data(n=100, er=0.6)
#'
#' # Split the data by 7/3
#'
#' subject_list <- unique(dat$Z$subjectid)
#'
#' train_id <- sample(subject_list, round(length(subject_list)*0.7))
#' test_id <- subject_list[!subject_list %in% train_id]
#'
#' # Data for training, keep only at risk data points
#' dat_train <- dat$Z %>%
#'     filter(at_risk == 1) %>%
#'     filter(subjectid %in% train_id) %>%
#'     select(-at_risk)
#'
#' fit <- dsp_fit(dat_train,
#'                C=2^2,
#'                kernel="semi_linear")
#'
#' # Data for testing
#' newdata <- dat$Z %>%
#'     filter(subjectid %in% test_id) %>%
#'     select(-at_risk, -event)
#'
#' # Making prediction
#' predict.dsp <- dsp_predict(fit, newdata = newdata)
#'
#'
#'
#'
#'
#' ## Example 2 - Make event time prediction for large testing data
#'
#' set.seed(1236)
#' dat <- sim_data(n=100, er=0.6)
#'
#' dat_train <- dat$Z %>%
#'     filter(at_risk == 1) %>%
#'     select(-at_risk)
#'
#' fit <- dsp_fit(dat_train,
#'                C=2^2,
#'                kernel="semi_linear")
#'
#'
#' dat_test <- sim_data(n=1000, er=1, time_points=fit$train_dat$unq_t_obs)
#'
#' predict.dsp <- dsp_predict(fit, select(dat_test$Z, -at_risk, -event))
#'




dsp_predict <- function(fit, newdata){

  test <- transform_Test(newdata, fit)


  train <- fit$train_dat

  if(any(class(fit) == "wsvm")){
    g <- predict(fit, train$x)
    g <- as.numeric(g) - 1
  }else{
    g <- predict_status(fit, x_test=train$x, x=train$x, kernel = fit$kernel, kparam = fit$sigma)
  }

  ## Get g downweight probability for training
  g_train_pred <-  data.frame(g = g,
                              at_risk_j_idx = train$at_risk_j_idx,
                              death = train$y) %>%
    filter(death == -1) %>%
    mutate(pred = ifelse(g > 0, 1, 0)) %>%
    group_by(at_risk_j_idx) %>%
    summarise(p_ = mean(pred))

  g_train_pred_event <-  data.frame(g = g,
                                    at_risk_j_idx = train$at_risk_j_idx,
                                    death = train$y) %>%
    filter(death == 1) %>%
    mutate(pred = ifelse(g > 0, 1, 0)) %>%
    group_by(at_risk_j_idx) %>%
    summarise(p_event = mean(pred))

  g_train_pred <- g_train_pred %>% left_join(g_train_pred_event, by="at_risk_j_idx")



  ## Get g for test data

  if(any(class(fit) == "wsvm")){
    g.test <- predict(fit, test$x)
    g.test <- as.numeric(g.test) - 1
  }else{
    g.test <- predict_status(fit, x_test = test$x, x=train$x, kernel = fit$kernel, kparam = fit$sigma)
  }

  n_test <- length(unique(test$subjectid))
  n_unq_t <- length(train$unq_t_obs)
  j_idx <- unlist(lapply(unique(test$subjectid),
                         function(x, i_idx){return(seq(length(i_idx[i_idx == x])))},
                         i_idx = test$subjectid
                         )
                  )

  test_pred <- data.frame(g=g.test,
                          subjectid=test$subjectid,
                          j_idx = j_idx,
                          event_time = train$unq_t_obs[j_idx],
                          sum = train$n_at_risk[j_idx],
                          n_death = train$n_death_t[j_idx]) %>%
    full_join(data.frame(p_ = g_train_pred$p_,
                         j_idx = g_train_pred$at_risk_j_idx,
                         p_event = g_train_pred$p_event), by=c("j_idx")) %>%
    mutate(p_ = ifelse(is.na(p_), 0, p_),
           p = 1/(p_*(sum/n_death-1)+1))



  rand_seq <- function(k, test_pred){
    test_pred_ <- test_pred
    test_pred_$r <- runif(dim(test_pred)[[1]])

    n_test <- length(unique(test_pred$subjectid))
    max_t <- max(test_pred$event_time)


    test_pred_ <- test_pred_ %>%
      mutate(event_ind = ifelse(r < p & g ==1, 1, 0)) %>%
      filter(event_ind == 1)

    test_pred_results <- test_pred_ %>%
      group_by(subjectid) %>%
      summarise(j_idx = min(j_idx)) %>%
      left_join(test_pred_, by=c("subjectid", "j_idx")) %>%
      right_join(data.frame(subjectid=unique(test_pred$subjectid)), by="subjectid") %>%
      mutate(event_time_final = ifelse(is.na(event_time), max_t, event_time),
             event_ind_final = ifelse(is.na(event_ind), 0, 1)) %>%
      arrange(subjectid) %>% dplyr::select(c(subjectid, event_time_final, event_ind_final))
    return(test_pred_results)
  }

  test_pred_results <- lapply(seq(50), rand_seq, test_pred=test_pred)
  test_pred_results_ <- do.call("bind_rows", test_pred_results)
  test_pred_results_ <- test_pred_results_ %>% group_by(subjectid) %>%
    summarise(event_time_final = mean(event_time_final))

  return(test_pred_results_)
}


predict_status <- function(fit, x_test, x, kernel, kparam){
  if(kernel == "semi_gaussian"){
    g_train <- apply(x_test,
                     MARGIN=1,
                     FUN=function(test, train, kernel, kparam, fit){
                            xinner <- xinner_kernel(x = test, y = train, kernel = kernel, kparam = kparam)
                            return(xinner %*% fit$beta + fit$beta0)
                            },
                     train=x,
                     kernel=kernel,
                     kparam=kparam,
                     fit=fit)

    g_status <- ifelse(g_train > 0, 1, 0)
  }else if(kernel == "semi_linear"){
    beta_train <- as.vector(colSums(as.numeric(fit$beta) * x)[-dim(x)[[2]]])
    beta_train_t <- aggregate(fit$beta, by= list(time=x[, dim(x)[[2]]]), FUN=sum)
    g_train <- as.matrix(x_test[, -dim(x)[[2]]]) %*% beta_train + beta_train_t$V1[match(x_test[, dim(x)[[2]]], beta_train_t$time)] + fit$beta0
    g_status <- ifelse(g_train > 0, 1, 0)
  }


  return(g_status)
}

