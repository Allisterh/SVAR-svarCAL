#install package
#' Demo model
#'
#' Basic RBC model to use as a testing function. The model works through the package gEcon make sure this package is installed before running the function.
#'
#' @param params list of parameters
#' @param T periods to simulate
#' @return
#' @export
#'
#' @examples
demo_model <- function(T, params){

  # make and load the model
  rbc <- make_model("R/rbc.gcn")

  delta <- params[1]
  beta <- params[2]
  eta <- params[3]
  mu <- params[4]
  phi <- params[5]
  rbc <-suppressWarnings(set_free_par(rbc, list(eta = eta, mu = mu, beta = beta, delta = delta, phi = phi)))

# find and print steady-state values
rbc <- steady_state(rbc)
get_ss_values(rbc, to_tex = F)

#find and print perturbation solution
rbc <- solve_pert(model = rbc, loglin = TRUE)
get_pert_solution(rbc, to_tex = F)

rbc <- set_shock_cov_mat(rbc, cov_matrix = matrix(0.01, 1, 1),
                         shock_order = "epsilon_Z")
rbc_ic_sim <- random_path(rbc, variables = c("K_s", "C", "Z", "I", "Y"),
                             sim_length = T)
return(rbc_ic_sim@sim)
}


