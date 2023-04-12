#' Model data generating function
#'
#' Function to generate the model data for multiple samples of the parameter space and multiple runs.
#'
#' @param model Function representing the theoretical model to be calibrated.
#' @param param_space A matrix consisting of two columns with the minimum and maximum value for each parameter.
#' @param sample_number Number of samples to draw from the parameter space.
#' @param runs The number of runs per parameter sample.
#' @param T Number of time periods to generate.
#'
#' @return A list of list of data frames containing the model data for the specified number of runs and samples.
#' @export
model_data <- function(model, param_space, sample_number, runs, T){
  # Make sure 'bounds' is a matrix with 2 columns
  if (ncol(param_space) != 2) {
    stop("param_space must have 2 columns")
  }

  # Generate 'samples' number of random samples from the sample spaces defined by the param space
  samples <- apply(param_space, 1, function(x) {
    runif(sample_number, min = x[1], max = x[2])
  })


  data_list <- list()
  data_list <- apply(samples, 1, function(y){
    tryCatch(
      lapply(runs, function(z){
        model(T,y)}), error = function(err) NA)

  })

  return(data_list)
}




#' Sample size function
#'
#' Function to calculate the minimum sample size following Seri & Secchi (2017).
#'
#' @param sample_number Number of parameter samples used.
#' @param effect_size Effect size, set to conservative default value of 0.1.
#' @param sigma significance level, set to default of 0.01
#' @param power_test power of ANAVO test, set to default of 0.95
#'
#' @return Numerical value indicating the minimum sample size given the inputs.
#'
#' @examples
#' sample_size(100)
#'
#' @export
sample_size <- function(sample_number , effect_size = 0.1, sigma = 0.01, power_test = 0.95){

  result <- pwr.anova.test(k = sample_number, n = NULL, f = effect_size, sig.level = sigma, power = power_test)
  return(ceiling(result$n))
}


