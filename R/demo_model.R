# ############################################################################
# (c) Chancellery of the Prime Minister 2012-2015                            #
#                                                                            #
# Authors: Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak         #
# ############################################################################
# Basic RBC model
# ############################################################################

# load gEcon package
library(gEcon)

# make and load the model
rbc <- make_model("R/rbc.gcn")

demo_model <- function(T, delta, beta, eta, mu, phi){

  rbc <- set_free_par(rbc, list(eta = eta, mu = mu, beta = beta, delta = delta, phi = phi))

# find and print steady-state values
rbc <- steady_state(rbc)
get_ss_values(rbc, to_tex = F)

# find and print perturbation solution
rbc <- solve_pert(model = rbc, loglin = TRUE)
get_pert_solution(rbc, to_tex = F)



rbc_ic_sim <- simulate_model(rbc, variables = c("K_s", "C", "Z", "I", "Y"),
                             shocks = c("epsilon_Z"),
                             shock_path = matrix(c(-0.05, 0, 0, -0.05), nrow = 1, ncol = 4),
                             sim_length = T)
return(rbc_ic_sim@sim)
}
test <- demo_model(50, 0.5, 0.5,0.5,0.5,0.5)
