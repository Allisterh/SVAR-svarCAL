#' Residual function of model estimated VAR
#'
#' Calculate the residuals based on the VAR estimate of the model data.
#'
#' @param var_est List of list of VAR estimates based on model data.
#'
#' @return List of lists of dataframes of residuals of VAR estimation
#' @export
resid_model <- function(var_est) {
  var_est <- m_VAR
  residuals <-lapply(var_est, function(x){
    lapply(x, function(y){
      q <- residuals(y)
     lapply(1:ncol(q), function(p){
       if(sum(q[,p]) == 0){
        q[1,p] <- 0.01}
     })
     q
    })

  })
  return(residuals)
}

#' Model fastICA function
#'
#' Identification of the list of list model mixing matrix through Indepdent component analysis alogrithm "fastICA".
#'
#' @param residuals List of lists of dataframes of VAR residuals
#' @inheritDotParams fastICA::fastICA -X -n.comp

#'
#' @return list of list of identified model mixing matrices
#' @export
fastICA_model <- function(residuals, ...){
  mix_matrix <-lapply(residuals, function(x){
    lapply(x, function(y){
      tryCatch(
      fastICA_gen(y,...) , error = function(err) NA)
    })
  })
  return(mix_matrix)
}


#' General fastICA function
#'
#' Indentification of the mixing matrix through Indepdent component analysis alogrithm "fastICA" and ordering of matrix.
#'
#' @param ures residuals of VAR model
#' @inheritDotParams fastICA::fastICA -X -n.comp
#'
#' @return Identified model mixing matrix.
#' @export
fastICA_gen <-function(ures, ...){
  set.seed(sseed)
  n<-ncol(ures)
  icares <- fastICA(ures, ncol(ures), ...)
  A<- t(icares$A)


  aba_select <- abs(A)
  n <- ncol(ures)
  a_new <- matrix(NA, nrow = n, ncol = n)
  for(i in 1:n){
    max_cor  <- which(aba_select == max(aba_select), arr.ind = TRUE)
    if(A[max_cor[1],max_cor[2]] < 0){
      a_new[, max_cor[1]] <- -A[,max_cor[2]]
    }
    else{
      a_new[, max_cor[1]] <- A[,max_cor[2]]
    }
    aba_select[, max_cor[2]] <- 0
    aba_select[max_cor[1],] <- 0
  }
  cnames<-paste("s", 1:n, sep="")
  return(a_new)
}



#' Moving average components function
#'
#'Function to extract the coefficients of lagged components of the moving average representation of the VAR model.
#'
#' @param var_est Result of estimated VAR model.
#' @param Mixmat Estimated mixing matrix of contemporaneous effects.
#' @param maxlag Number of lagged coefficients to compute, default is set to 5.
#'
#' @return List of laged component matrices.
#' @export
lag_MA_components <- function(var_est, Mixmat, maxlag = 5){
  AA <- Acoef(var_est)
  k <- nrow(AA[[1]])
  p <- length(AA)
  AAs <- as.list(1:(maxlag))
  AA <- Acoef(var_est)
  AAs <-lapply(1:(maxlag), function(x){
    matrix(0,nrow=k,ncol=k)
  })
  AAs[1:p] <- AA

  Fi <- lapply(1:(maxlag+1), function(i) { if (i == 1) diag(k) else if (i == 2) AAs[[1]] else matrix(0, nrow = k, ncol = k)})
  for (i in 2:maxlag){
    for (j in 1:i){
      FF<- Fi[[i+1-j]]%*%AAs[[j]]
      Fi[[i+1]]<-Fi[[i+1]]+FF
    }
  }

 lapply(Fi, function(x){x <- x %*% Mixmat})
}

########


#' lag_MA_components_mod
#'
#'Function to apply the lag_MA_components function to a list of list of model estimates
#' @param var_est_model List of list of VAR estimates of the model.
#' @param Mixmat_model List of lists of estimates of the model mixing matrix.
#' @inheritDotParams lag_MA_components maxlag
#'
#' @return List of lists of lists of lagged component matrices.
#' @export
lag_MA_components_mod <- function(var_est_model, Mixmat_model, ...){
  p <- (var_est_model[[1]][[1]][['p']]+1)
  result <- Map(function(a, b) {
    Map(function(c, d) {
      if (!any(is.na(d))) {
        lag_MA_components(c, d, ...)
      } else {
        rep(NA,p)
      }
    }, a, b)
  }, var_est_model, Mixmat_model)
  return(result)
}



