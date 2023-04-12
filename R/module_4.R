#' Minimal distance index function
#'
#'Function to calculate the Minimal Dinstance Index (MDI) between the estimates of the model mixing matrix and the real world mixing matrix. The distance measure is the frobenius norm.
#'
#' @param phi_mod MA representation of model data VAR
#' @param phi_emp MA representation of empirical data VAR
#'
#' @return List of lagged components each consisting of a list with parameter settings consisting of dataframes with frobenius distances for each run.
#' @export
MDI <- function(phi_mod, phi_emp) {

  K <- ncol(phi_emp[[1]])
  L <- length(phi_emp)


  #calculate MDI
  frob_value <- new.env()
  Q <- lapply(phi_mod, function(q) {

    lapply(q, function(x) {


      frob_value$t <- list()

      lapply(1:length(x), function(y) {

        if (length(x[[y]]) == K*K & length(phi_emp) == length(x) & !any(is.na(x[[y]]))) {
          if (y == 1) {

              frob_value$t <- signperm_matrix(x[[y]], phi_emp[[y]]) #compute optimal sign perm matrix
              frobenius.norm(x[[y]] %*% frob_value$t - phi_emp[[y]]) #what about (1/sqrt(K-1))

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


#' Sign permutation function
#'
#' Compute optimal sign permutation matrix using linear assignment algorithm.
#'
#' @param A n by n matrix such that ||AP-B||f.
#' @param B n by n matrix such that ||AP-B||f.
#'
#' @return Sign permutation matrix that minimizes the frobenius norm between A and B
#' @export
signperm_matrix <- function(A, B) {
  D_plus <- t(A) %*% B #compute cost matrix D if A is positive, each matrix entry is the multiplication between a column of A and a column of B (one entry for each combination of A and B)
  D_min <- t(-A) %*% B #compute cost matrix D  if A is negative
  D <- D_plus*(D_plus >= D_min) + D_min * (D_min > D_plus) # select those elements of D_plus and D_min that are the highest into a new matrix D
  sign_change <- 1 * (D_plus >= D_min) + (-1) * (D_min > D_plus) #track which elements are from D_plus and D_minus
  P <- lp.assign(-D)$solution #use linear assignment algo to determine permutation matrix that maximizes the cost matrix (algorithm minimizes by default so -D)
  perm_sign <- P *  sign_change #change sign of permutation matrix if element came from D_min
  return(perm_sign)
}


#' P values function
#'
#' Function to determine the p-values of minimal distance index based on Chi-squared. This function is used within the MCS() function
#'
#' @param iN Number of runs per CoP.
#' @param dA Significance theshold.
#' @param D_bar Vector of average distances for different CoPs.
#' @param variance Covariance matrix of D_bar.
#'
#' @return List of results about the model confidence set
#' @export
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
#' @return Model confidence set information
#' @export
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
#' @param fun_el elimination rule
#'
#' @return Model confidence set of one lagged component.
#' @export
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
#' @inheritDotParams MCS -MDI_matrix
#'
#' @return List of lagged components consisting of the model confidence sets.
#' @export
MCS_list <- function(MDI_matrix_list, ...){

  lapply(MDI_matrix_list, function(x){
    MCS(dA = dA, MDI_matrix = x, fun_p = p_values, fun_el = elim)
  })
}


