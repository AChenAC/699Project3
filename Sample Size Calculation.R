###Simulations###
trial_sim_bernoulli <- function(n_sim, 
                                n_A, 
                                n_B, 
                                pA,
                                pB,
                                seed = sample(.Machine$integer.max, 1)) {
  set.seed(seed)
  start_sim = Sys.time();
  all_data_A <- 
    rbinom(n_sim, n_A, pA) / n_A;
  all_data_B <-
    rbinom(n_sim, n_B, pB) / n_B;
  estimated_sigmasq_A <- all_data_A * (1 -  all_data_A);
  estimated_sigmasq_B <- all_data_B * (1 -  all_data_B);
  t_stat <-
    (all_data_A - all_data_B) / 
    sqrt(estimated_sigmasq_A / n_A + estimated_sigmasq_B / n_B);
  df_approx <- 
    (estimated_sigmasq_A/n_A + estimated_sigmasq_B/n_B)^2 / 
    ((estimated_sigmasq_A/n_A)^2/(n_A-1) + 
       (estimated_sigmasq_B/n_B)^2/(n_B-1));
  list(pvalues_t = 2 * pt(q = abs(t_stat), 
                          df = df_approx, 
                          lower.tail = FALSE), 
       runtime = Sys.time() - start_sim);  
}



#Interim analysis: Under the null
pvals_interim_null_1 = unname(unlist(trial_sim_bernoulli(2000, 59, 59, 0.15, 0.15, seed = 456)))
length(which(pvals_interim_null_1<0.005))/length(pvals_interim_null_1)

pvals_interim_null_2 = unname(unlist(trial_sim_bernoulli(2000, 59, 59, 0.05, 0.05, seed = 456)))
length(which(pvals_interim_null_2<0.005))/length(pvals_interim_null_2)

#Type 1 error proportion is <0.05 for each scenario.

#Interim analysis: Under the alternative
pvals_interim_alt_1 = unname(unlist(trial_sim_bernoulli(2000, 59, 59, 0.15, 0.3, seed = 456)))
length(which(pvals_interim_alt_1<0.005))/length(pvals_interim_alt_1)

pvals_interim_alt_2 = unname(unlist(trial_sim_bernoulli(2000, 59, 59, 0.05, 0.2, seed = 456)))
length(which(pvals_interim_alt_2<0.005))/length(pvals_interim_alt_2)

#Power is 0.20 for scenario 1 and 0.36 for scenario 2. This means we are unlikely to
#detect a true difference in the interim analysis with this sample size.



#Final analysis: Under the null

pvals_final_null_1 = unname(unlist(trial_sim_bernoulli(2000, 119, 119, 0.15, 0.15, seed = 456)))
length(which(pvals_final_null_1<0.048))/length(pvals_final_null_1)


pvals_final_null_2 = unname(unlist(trial_sim_bernoulli(2000, 119, 119, 0.05, 0.05, seed = 456)))
length(which(pvals_final_null_2<0.048))/length(pvals_final_null_2)

#Final analysis: Under the alternative
pvals_final_alt_1 = unname(unlist(trial_sim_bernoulli(2000, 119, 119, 0.15, 0.3, seed = 456)))
length(which(pvals_final_alt_1<0.048))/length(pvals_final_alt_1)


pvals_final_alt_2 = unname(unlist(trial_sim_bernoulli(2000, 119, 119, 0.05, 0.2, seed = 456)))
length(which(pvals_final_alt_2<0.048))/length(pvals_final_alt_2)
