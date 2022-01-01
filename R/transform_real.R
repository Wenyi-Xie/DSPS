data_transform <- function(base_data){

  base_data <- as.data.frame(base_data)
  # Event occurrence
  first_occ <- base_data %>% filter(event == 1)

  # Unique event time
  unq_t_obs <- sort(unique(first_occ$time))

  temp <- base_data %>% filter(time %in% c(unq_t_obs))


  # Subject level last visit status
  base_data_ <- temp %>% group_by(subjectid) %>% summarise(T_obs=max(time),event_ind_obs = max(event))

  base_data_ <- left_join(data.frame(subjectid = unique(temp$subjectid)), base_data_, by="subjectid")
  expanded_Z <- temp %>% dplyr::select(-event, -subjectid)


  T_obs_train <- base_data_$T_obs

  event_ind_obs_train <- base_data_$event_ind_obs





  ### ---- Data preparation ----
  n <- length(unique(base_data_$subjectid)) # Total training sample size n

  # Total number of unique event time
  n_unq_t <- length(unq_t_obs)

  # Event time index for subjects failed
  At_risk_time_length <- sapply(T_obs_train, FUN=function(x){sum(x>=unq_t_obs)})

  # At risk process
  Y_at_risk_t <- matrix(0, nrow=n_unq_t, ncol=n)
  idx_incidence <- ifelse(event_ind_obs_train==1, match(T_obs_train, unq_t_obs), NA)

  Y_at_risk_t[unlist(
    sapply(1:n, function(x){
      if(At_risk_time_length[x]==0){return(NA)}
      else{return(((x - 1) * n_unq_t + 1):((x - 1) * n_unq_t + At_risk_time_length[x]))}
    }
    ))] <- 1

  # Y_at_risk is of n by n_unq_t
  Y_at_risk <- t(Y_at_risk_t)

  ## complete_vector_at_risk is Y_at_risk expanded
  ## length: n * n_unq_t
  complete_vector_at_risk <- as.vector(t(Y_at_risk))

  # Counting process
  delta_N <- matrix(data=-1, nrow=n, ncol = n_unq_t)
  delta_N[cbind(1:n, idx_incidence)] <- 1 # delta_N of dim n by n_unq_t

  ## Expanded counting process with only at risk points
  at_risk_only_delta <- as.vector(t(delta_N))[complete_vector_at_risk==1]


  # i and j tracking index
  at_risk_i_idx <- rep(1:n, each=n_unq_t)[complete_vector_at_risk==1]
  at_risk_j_idx <- rep(1:n_unq_t, n)[complete_vector_at_risk==1]

  n_at_risk_i_j <- length(at_risk_i_idx)



  # Weight w
  weight_fn <- function(i){
    weight_fn_1 <- function(i, j){
      death_tj <- sum(delta_N[, j]==1)
      if(delta_N[i,j] == 1){
        at_risk_total_tj <- sum(Y_at_risk[, j])
        return(1 - death_tj/at_risk_total_tj)
      }else if(Y_at_risk[i,j]==1){
        return(death_tj/sum(Y_at_risk[, j]))
      }#else{
      #return(NA)
      #}
    }
    return(unlist(sapply(1:n_unq_t, weight_fn_1, i=i)))
  }


  wi_tj <- unlist(sapply(1:n, weight_fn))

  return(list(weight = wi_tj,
              y = at_risk_only_delta,
              x = expanded_Z,
              at_risk_i_idx = at_risk_i_idx,
              at_risk_j_idx = at_risk_j_idx,
              unq_t_obs = unq_t_obs,
              n_at_risk = colSums(Y_at_risk),
              n_death_t = colSums(delta_N == 1),
              id = temp$subjectid))
}


transform_Test <- function(test, fit){

  temp <- test %>% filter(time %in% fit$train_dat$unq_t_obs)

  Z_test_expanded <- temp %>% dplyr::select(-subjectid)

  Z_test_expanded <- sweep(sweep(Z_test_expanded, 2, fit$train_mean), 2, fit$train_sd, FUN="/")



  subjectid <- temp$subjectid

  return(list(x=Z_test_expanded, subjectid = subjectid))
}

