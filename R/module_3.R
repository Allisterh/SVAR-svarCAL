#' resid_model
#'
#' Calculate the residuals based on the VAR estimate of the model data
#'
#' @param var_est VAR estimate based on model data
#'
#' @return
#' @export
#'
#' @examples
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

#' fast_ICA_model
#'
#' Indentification of the list of list model mixing matrix through Indepdent component analysis alogrithm "fastICA"
#'
#' @param residuals matrix of VAR residuals
#' @param sseed seed
#' @param tol see fastICA docu
#' @param maxit see fastICA docu
#' @param verbose see fastICA docu
#'
#' @return
#' @export
#'
#' @examples
fastICA_model <- function(residuals, sseed = 46, maxit=3000, tol=1e-14 , verbose=FALSE){
  mix_matrix <-lapply(residuals, function(x){
    sseed=46
    lapply(x, function(y){
      tryCatch(
      fastICA_gen(y, sseed=46) , error = function(err) NA)
    })
  })
  return(mix_matrix)
}


#' fast_ICA_general
#'
#' Indentification of the mixing matrix through Indepdent component analysis alogrithm "fastICA" and ordering of matrix
#' @param ures residuals of VAR model
#' @param sseed seed
#' @param tol see fastICA docu
#' @param maxit see fastICA docu
#' @param verbose see fastICA docu
#'
#' @return
#' @export
#'
#' @examples
fastICA_gen <-function(ures, sseed=46, maxit=3000, tol=1e-14 , verbose=FALSE){
  set.seed(sseed)
  n<-ncol(ures)
  icares <- fastICA(ures, ncol(ures),tol=tol, maxit=maxit, verbose=verbose) #Doris: tol=1e-14
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



#' lag_MA_components
#'
#'Function to extract the coefficients of  lagged componenets of the moving average represenetation of the VAR model
#'
#' @param var_est result of estimated VAR model
#' @param Mixmat estimated mixing matrix of contemporanous effects
#' @param maxlag nubmer of lagged coefficient to compute
#'
#' @return
#' @export
#'
#' @examples
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
#' @param var_est_model var estimates of the model in a list of lists
#' @param Mixmat_model estimates of the model mixing matrx in a list of lists
#' @param maxlag number of lagged components to compute
#'
#' @return
#' @export
#'
#' @examples
lag_MA_components_mod <- function(var_est_model, Mixmat_model, maxlag = 5){
  p <- (var_est_model[[1]][[1]][['p']]+1)
  result <- Map(function(a, b) {
    Map(function(c, d) {
      if (!any(is.na(d))) {
        lag_MA_components(c, d, maxlag = maxlag)
      } else {
        rep(NA,p)
      }
    }, a, b)
  }, var_est_model, Mixmat_model)
  return(result)
}



