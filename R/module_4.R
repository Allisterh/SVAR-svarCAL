
#space <- rbind(c(0.99*0.9, 0.99*1.1), c(0.9*0.025, 0.025*1.1), c(6.771*0.9, 6.771*1.1), c(0.169*1, 1*0.169), c(1*1.353, 1*1.353), c(1*0.573,1*0.573), c(1*2.4,1*2.4), c(1*1,1*1))
#test <- model_data(model = demo_model, param_space = space, sample_number = 5, runs = 10, T = 200)
#data("Philippine_SVAR")

#var_model_est <- VAR_model(model_data = test, var_select = c('Y', 'I', 'pi'), trim = 10, lag_select = 4, season = NULL, exog = NULL, type = 'const' )
#var_emps_est <- VAR_emp(data = Philippine_SVAR, var_select = c('Output.Gap', 'RRP', 'CPI'), lag_select =4, season = NULL, exog = NULL, type = 'const' )
#input <- resid_model(var_model_est)
#test_mod <- fast_ICA_model(input, n.comp = 3)
#test_emp <- fastICA(test_4, n.comp = 2)


