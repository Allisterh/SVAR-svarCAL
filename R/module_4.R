#' MDI
#'
#'Function to calculate the Minimal Dinstance Index (MDI) between the estimates of the model mixing matrix and the real world mixing matrix
#' @param phi_mod MA representation of model data VAR
#' @param phi_emp MA representation of empirical data VAR
#'
#' @return
#' @export
#'
#' @examples
MDI <- function(phi_mod, phi_emp) {

  K <- ncol(phi_emp[[1]])
  L <- length(phi_emp)


  #create set of all P
  g <- diag(K)
  perms <- permutations(K,K)
  pairings <- expand.grid(1:nrow(perms))
  P_all <- lapply(1:nrow(pairings), function(x) g[perms[pairings[x,1],],])

  #create set of all J
  signs <- expand.grid(replicate(K, c(1,-1), simplify=FALSE))
  J_all <- lapply(1:nrow(signs), function(x) matrix(diag(signs[x,]), ncol = K))


  #create set of all C
  C_all <- flatten(lapply(P_all, function(x) lapply(J_all, function(y) x %*% y )))


  #calculate MDI
  frob_value <- new.env()
  Q <- lapply(phi_mod, function(q) {

    lapply(q, function(x) {


      frob_value$t <- list()

      lapply(1:length(x), function(y) {

        if (length(x[[y]]) == K*K & length(phi_emp) == length(x) & !any(is.na(x[[y]]))) {
          if (y == 1) {
            result <- sapply(C_all, function(z) {

              frobenius.norm(x[[y]] %*% z - phi_emp[[y]]) #what about (1/sqrt(K-1))?

            })
            frob_value$t <- C_all[[which.min(result)]] #store sign order matrix that minimizes result

            min(result)

          } else {
            frobenius.norm(x[[y]] %*% frob_value$t - phi_emp[[y]]) #caluclate MDI for the lags
          }

        } else {
          NA
        }
      })
    })
  })

    #Transpose lists such that we get a list of one dataframe for each lag
    Q <- lapply(Q, function(x){
        purrr::transpose(x)
    })
    Q <- purrr::transpose(Q)
    names(Q) <- paste('phi', seq(from = 0, to = (L-1), by = 1))

    Q <- lapply(Q, function(x){
      lapply(1:length(x), function(y){
        x[[y]] <- as.data.frame(unlist(x[[y]]))
        names(x[[y]]) <- names(x)[[y]]
        x[[y]]
      })
    })
    Q <- lapply(Q, function(x){
      do.call(cbind, x)
    })

  return(Q)
}



#' p_values
#'
#' Function to determine the p-values of minimal distance index based on Chi-squared
#'
#'
#' @param iN number of runs per CoP
#' @param dA significance theshold
#' @param D_bar vector of average distances for different CoPs
#' @param variance covariance matrix of D_bar
#'
#' @return
#' @export
#'
#' @examples
p_values <- function(iN, dA, D_bar, variance){
  #D_bar
  iM <- length(D_bar)


  if(iM>2){
    #sigma and sigma_1
    Sigma <- diag(variance)
    Sigma_1 <- Sigma[1,1]
    Sigma <- as.matrix(Sigma[-1,-1])
    Sigma_n <- Sigma + Sigma_1

    #matrix
    A <- cbind(rep(1,iM -1),diag(-rep(1,iM -1)))
    test <- iN*t(A%*%D_bar) %*% solve(Sigma_n) %*% (A%*%D_bar)

  } else {
    test <- iN*(D_bar[1] -D_bar[2])**2/(variance[1] +variance[2])
  }

  if(iM>1){
    dQuan <- qchisq(1-dA,iM-1)
    dPVa <- 1-pchisq(test,iM-1)
  } else {
    dQuan <- 0
    dPVa <- 1
  }
  dResult <- 1*(test>dQuan)
  if (is.na(dResult)) dResult <- 0
  return(list(test=test,q=dQuan,p=dPVa,result=dResult, D_bar = D_bar, variance = variance))
}



#elimination function
#' elim
#'
#' @param D_bar vector of average distances for different CoPs
#' @param variance covariance matrix of D_bar
#'
#' @return
#' @export
#'
#' @examples
elim <- function(D_bar, variance){
  ind <- which.max(D_bar)
  MDI_max <- max(D_bar)
  variance_max <- variance[names(ind)]
  return(list(variance = variance[names(variance) != names(ind)],D_bar = D_bar[names(D_bar) != names(ind)], name = names(ind), MDI = MDI_max, var = variance_max))
}



#MDI assessment function
#' MCS function
#'
#' Function to determine the model confidence set and p-values
#'
#' @param dA significance threshold
#' @param MDI_matrix matrix of MDIs
#' @param fun_p function to determine p-value
#' @param fun_el elemination rule
#'
#' @return
#' @export
#'
#' @examples
MCS <- function(MDI_matrix, fun_p = p_values, fun_el = elim, dA = dA){

  MDI_matrix <- t(MDI_matrix)
  n_run <- ncol(MDI_matrix)
  n_cop <- nrow(MDI_matrix)

  #determine initial D_bar
  D_bar <-rowMeans(MDI_matrix, na.rm = F)
  D_bar[is.na(D_bar)] <- 10000


  #determine initial variance
  variance <- rowVars(MDI_matrix, useNames = T)
  variance[is.na(variance)] <- 10000

  #initialize results matrix
  result_matrix <- as.data.frame(matrix(ncol = 6, nrow = n_cop))
  names(result_matrix) <- c('test', 'q', 'p', 'result', 'MDI', 'var')
  for(i in 1:(n_cop)){

    #run results function
    result <- fun_p(dA = dA, iN = n_run, D_bar = D_bar, variance = variance)

    #add to results matrix
    result_matrix[i,1:4] <- unlist(result[c('test', 'q', 'p', 'result')])

    #eliminate highest D_bar
    result_elim <- fun_el(result$D_bar, result$variance)
    row.names(result_matrix)[i] <- result_elim$name
    result_matrix[i,5] <- result_elim$MDI
    result_matrix[i,6] <- result_elim$var
    D_bar <- unlist(result_elim$D_bar)
    variance <- unlist(result_elim$variance)

  }
  return(result_matrix)
}


#' MSC_list
#'
#' Computes the minimal confidence set of list of MDI matrices
#'
#' @param MDI_matrix_list list of MDI matrices
#' @param dA significance level default dA=1 selet best cop
#'
#' @return
#' @export
#'
#' @examples
MCS_list <- function(MDI_matrix_list, dA = 1){

  lapply(MDI_matrix_list, function(x){
    MCS(dA = dA, MDI_matrix = x, fun_p = p_values, fun_el = elim)
  })
}


