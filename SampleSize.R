## Scenario 1
target_size <- 0.05;
target_power <- 0.80;
assumed_pA = 0.15;
assumed_pB = 0.30;
assumed_delta <- assumed_pB - assumed_pA;
assumed_sigmasq_A <- assumed_pA * (1 - assumed_pA);
assumed_sigmasq_B <- assumed_pB * (1 - assumed_pB);
ratio_n_B_n_A <- 1;
required_n_A <- (qnorm(1 - target_size/2) + qnorm(target_power))^2 / assumed_delta^2 * 
  (assumed_sigmasq_A + assumed_sigmasq_B / ratio_n_B_n_A);
required_n_B <- ratio_n_B_n_A * required_n_A;
# bigger assumed_pa needs larger sample size 

ceiling(119/2) # 60

# Interim
sim_delta_eq_0 = trial_sim_bernoulli(2000, 60, 60, 0.15, 0.15, 1) # Null
mean(sim_delta_eq_0$pvalues_t < 0.005);
sim_delta_eq_target = trial_sim_bernoulli(2000, 60, 60, 0.15, 0.30, 1) # Alternative
mean(sim_delta_eq_target$pvalues_t < 0.005);

# Final
sim_delta_eq_0 = trial_sim_bernoulli(2000, 119, 119, 0.15, 0.15, 1) # Null
mean(sim_delta_eq_0$pvalues_t < 0.048);
sim_delta_eq_target = trial_sim_bernoulli(2000, 119, 119, 0.15, 0.30, 1) # Alternative
mean(sim_delta_eq_target$pvalues_t < 0.048);

## Scenario 2
target_size <- 0.05;
target_power <- 0.80;
assumed_pA = 0.05;
assumed_pB = 0.20;
assumed_delta <- assumed_pB - assumed_pA;
assumed_sigmasq_A <- assumed_pA * (1 - assumed_pA);
assumed_sigmasq_B <- assumed_pB * (1 - assumed_pB);
ratio_n_B_n_A <- 1;
required_n_A <- (qnorm(1 - target_size/2) + qnorm(target_power))^2 / assumed_delta^2 * 
  (assumed_sigmasq_A + assumed_sigmasq_B / ratio_n_B_n_A);
required_n_B <- ratio_n_B_n_A * required_n_A;
# bigger assumed_pa needs larger sample size 

ceiling(73/2) # 37

# Interim
sim_delta_eq_0 = trial_sim_bernoulli(2000, 37, 37, 0.05, 0.05, 1) # Null
mean(sim_delta_eq_0$pvalues_t < 0.005);
sim_delta_eq_target = trial_sim_bernoulli(2000, 37, 37, 0.05, 0.20, 1) # Alternative
mean(sim_delta_eq_target$pvalues_t < 0.005);

# Final
sim_delta_eq_0 = trial_sim_bernoulli(2000, 37, 37, 0.05, 0.05, 1) # Null
mean(sim_delta_eq_0$pvalues_t < 0.048);
sim_delta_eq_target = trial_sim_bernoulli(2000, 37, 37, 0.05, 0.20, 1) # Alternative
mean(sim_delta_eq_target$pvalues_t < 0.048);

