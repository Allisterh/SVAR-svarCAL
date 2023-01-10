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
  sw_gecon <- make_model("R/SW_03.gcn")

  beta <- params[1] #0.99 discount factor
  tau <- params[2] # 0.025 Capital depreciation rate
  varphi <- params[3] # 6.771 Parameter of investment adjustment cost function
  psi <- params[4] # 0.169 Capacity utilisation cost parameter
  sigma_c <- params[5] #1.353 Coefficient of relative risk aversion
  h <- params[6] # 0.573 Habit formation intensity
  sigma_l <- params[7] # 2.4 Reciprocal of labour elasticity w.r.t. wage
  omega <- params[8] #1; Labour disutility parameter

  sw_gecon <- suppressWarnings(set_free_par(sw_gecon, list(beta = beta, tau = tau, varphi = varphi, psi = psi, sigma_c = sigma_c, h = h, sigma_l = sigma_l, omega = omega)))

  # set initial values
  initv <- list(z = 1, z_f = 1, Q = 1, Q_f = 1, pi = 1,
                pi_obj = 1, epsilon_b = 1, epsilon_L = 1,
                epsilon_I = 1, epsilon_a = 1, epsilon_G = 1,
                r_k = 0.01, r_k_f = 0.01)
  sw_gecon <- initval_var(sw_gecon, init_var = initv)

  # find and print steady-state values
  sw_gecon <- steady_state(sw_gecon)
  get_ss_values(sw_gecon, to_tex = TRUE)

  # find and print perturbation solution
  sw_gecon <- solve_pert(sw_gecon, loglin = TRUE)
  get_pert_solution(sw_gecon, to_tex = TRUE)

  # set the shock distribution parameters
  variances <- c(eta_b = 0.336 ^ 2, eta_L = 3.52 ^ 2, eta_I = 0.085 ^ 2,
                 eta_a = 0.598 ^ 2, eta_w = 0.6853261 ^ 2, eta_p = 0.7896512 ^ 2,
                 eta_G = 0.325 ^ 2, eta_R = 0.081 ^ 2, eta_pi = 0.017 ^ 2)
  sw_gecon  <- set_shock_cov_mat(sw_gecon,
                                 cov_matrix = diag(variances),
                                 shock_order = names(variances))


rbc_ic_sim <- random_path(sw_gecon, variables = c("q", "pi", "r_k",
                                                  "z", "C", "G", "I", "K", "L",
                                                  "Q", "R", "W", "T", "Y"),
                             sim_length = T)
return(rbc_ic_sim@sim)
}

