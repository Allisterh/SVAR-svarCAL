#' Model Data Function
#'
#' Function to generate the model data for multiple samples of the parameter space and multiple runs
#'
#' @param model function representing the theoretical model to be calibrated
#' @param param_space a matrix consisting of two columns with the minimum and maximum value for each parameter
#' @param sample_number number of samples to draw from the parameter space
#' @param runs the number of runs per parameter sample
#' @param T number of time periods to generate
#'
#' @return a matrix containing the model data for the specified number of runs and samples
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
      lapply(1:20, function(z){
        model(T,y)}), error = function(err) NA)

  })

  return(data_list)
}

#space <- rbind(c(0.025*0.99, 0.025*1.01), c(0.99*0.99, 0.99*1.01), c(2*0.99, 2*1.01), c(0.3*1, 1*0.3), c(1*0.95, 1*0.95))
#test <- model_data(model = demo_model, param_space = space, sample_number = 5, runs = 10, T = 50)


#' Sample size
#'
#' Function to calculate sample following Seri & Secchi (2017)
#'
#' @param sample_number number of parameter samples used
#' @param effect_size set to conservative default of 0.1
#' @param sigma significance level, set to default of 0.01
#' @param power_test power of ANAVO test, set to default of 0.95
#'
#' @return
#' @export
#'
#' @examples
sample_size <- function(sample_number , effect_size = 0.1, sigma = 0.01, power_test = 0.95){

  result <- pwr.anova.test(k = sample_number, n = NULL, f = effect_size, sig.level = sigma, power = power_test)
}


