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
  residuals <-lapply(var_est, function(x){
    lapply(x, function(y){
      residuals(y)

    })
  })
  return(residuals)
}

#' fast_ICA_model
#'
#' Indentification of the model mixing matrix through Indepdent component analysis alogrithm "fastICA"
#'
#' @param residuals matrix of VAR residuals
#' @param n.comp number of components to identify
#' @param ... see other fastICA comments
#'
#' @return
#' @export
#'
#' @examples
fast_ICA_model <- function(residuals, n.comp, ...){
  mix_matrix <-lapply(residuals, function(x){

    lapply(x, function(y){
      tryCatch(
      fastICA(y, n.comp) , error = function(err) NA)
    })
  })
  return(mix_matrix$K)
}




