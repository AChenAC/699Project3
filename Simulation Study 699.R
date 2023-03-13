library(tidyverse)


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



###Trying different sample sizes###

#Interim analysis
pvals_interim_null_1 = list()
pvals_interim_null_2 = list()
pvals_interim_null_3 = list()
type1error_interim_1 = c()
type1error_interim_2 = c()
type1error_interim_3= c()

pvals_interim_alt_1 = list()
pvals_interim_alt_2 = list()
pvals_interim_alt_3 = list()
power_interim_1 = c()
power_interim_2 = c()
power_interim_3 = c()

for(i in 110:130){
  
  pvals_interim_null_1[[i]] = unname(unlist(trial_sim_bernoulli(2000, i, i, 0.05, 0.05, seed = 456)))
  pvals_interim_null_2[[i]] = unname(unlist(trial_sim_bernoulli(2000, i, i, 0.1, 0.1, seed = 456)))
  pvals_interim_null_3[[i]] = unname(unlist(trial_sim_bernoulli(2000, i, i, 0.15, 0.15, seed = 456)))
  
  type1error_interim_1[i] = length(which(pvals_interim_null_1[[i]]<0.005))/length(pvals_interim_null_1[[i]])
  type1error_interim_2[i] = length(which(pvals_interim_null_2[[i]]<0.005))/length(pvals_interim_null_2[[i]])
  type1error_interim_3[i] = length(which(pvals_interim_null_3[[i]]<0.005))/length(pvals_interim_null_3[[i]])
  
  pvals_interim_alt_1[[i]] = unname(unlist(trial_sim_bernoulli(2000, i, i, 0.05, 0.2, seed = 456)))
  pvals_interim_alt_2[[i]] = unname(unlist(trial_sim_bernoulli(2000, i, i, 0.1, 0.25, seed = 456)))
  pvals_interim_alt_3[[i]] = unname(unlist(trial_sim_bernoulli(2000, i, i, 0.15, 0.3, seed = 456)))
  
  power_interim_1[i] = length(which(pvals_interim_alt_1[[i]]<0.005))/length(pvals_interim_alt_1[[i]])
  power_interim_2[i] = length(which(pvals_interim_alt_2[[i]]<0.005))/length(pvals_interim_alt_2[[i]])
  power_interim_3[i] = length(which(pvals_interim_alt_3[[i]]<0.005))/length(pvals_interim_alt_3[[i]])
}

type1error_interim_1 = type1error_interim_1[110:130]
type1error_interim_2 = type1error_interim_2[110:130]
type1error_interim_3 = type1error_interim_3[110:130]

power_interim_1 = power_interim_1[110:130]
power_interim_2 = power_interim_2[110:130]
power_interim_3 = power_interim_3[110:130]


#Final Analysis
pvals_final_null_1 = list()
pvals_final_null_2 = list()
pvals_final_null_3 = list()
type1error_final_1 = c()
type1error_final_2 = c()
type1error_final_3= c()

pvals_final_alt_1 = list()
pvals_final_alt_2 = list()
pvals_final_alt_3 = list()
power_final_1 = c()
power_final_2 = c()
power_final_3 = c()

for(i in 110:130){
  
  pvals_final_null_1[[i]] = unname(unlist(trial_sim_bernoulli(2000, i, i, 0.05, 0.05, seed = 456)))
  pvals_final_null_2[[i]] = unname(unlist(trial_sim_bernoulli(2000, i, i, 0.1, 0.1, seed = 456)))
  pvals_final_null_3[[i]] = unname(unlist(trial_sim_bernoulli(2000, i, i, 0.15, 0.15, seed = 456)))
  
  type1error_final_1[i] = length(which(pvals_final_null_1[[i]]<0.048))/length(pvals_final_null_1[[i]])
  type1error_final_2[i] = length(which(pvals_final_null_2[[i]]<0.048))/length(pvals_final_null_2[[i]])
  type1error_final_3[i] = length(which(pvals_final_null_3[[i]]<0.048))/length(pvals_final_null_3[[i]])
  
  pvals_final_alt_1[[i]] = unname(unlist(trial_sim_bernoulli(2000, i, i, 0.05, 0.2, seed = 456)))
  pvals_final_alt_2[[i]] = unname(unlist(trial_sim_bernoulli(2000, i, i, 0.1, 0.25, seed = 456)))
  pvals_final_alt_3[[i]] = unname(unlist(trial_sim_bernoulli(2000, i, i, 0.15, 0.3, seed = 456)))
  
  power_final_1[i] = length(which(pvals_final_alt_1[[i]]<0.048))/length(pvals_final_alt_1[[i]])
  power_final_2[i] = length(which(pvals_final_alt_2[[i]]<0.048))/length(pvals_final_alt_2[[i]])
  power_final_3[i] = length(which(pvals_final_alt_3[[i]]<0.048))/length(pvals_final_alt_3[[i]])
}

type1error_final_1 = type1error_final_1[110:130]
type1error_final_2 = type1error_final_2[110:130]
type1error_final_3 = type1error_final_3[110:130]

power_final_1 = power_final_1[110:130]
power_final_2 = power_final_2[110:130]
power_final_3 = power_final_3[110:130]

type1error_interim = data.frame(cbind(type1error_interim_1, type1error_interim_2, type1error_interim_3))
type1error_final = data.frame(cbind(type1error_final_1, type1error_final_2, type1error_final_3))

power_interim = data.frame(cbind(power_interim_1, power_interim_2, power_interim_3))
power_final = data.frame(cbind(power_final_1, power_final_2, power_final_3))

data = cbind(type1error_final, power_final, type1error_interim, power_interim)
final = cbind(type1error_final, power_final)
interim = cbind(type1error_interim, power_interim)

final$n = 110:130
interim$n = 110:130
data$n = 110:130

newfinal = final %>%
  pivot_longer(c(type1error_final_1, type1error_final_2, type1error_final_3, power_final_1, power_final_2, power_final_3),
               names_to = "effect_sizes", values_to = "value")

newinterim = interim %>%
  pivot_longer(c(type1error_interim_1, type1error_interim_2, type1error_interim_3, power_interim_1, power_interim_2, power_interim_3),
               names_to = "effect_sizes", values_to = "value")
newdata = data %>%
  pivot_longer(c(type1error_final_1, type1error_final_2, type1error_final_3, power_final_1, power_final_2, power_final_3, 
                 type1error_interim_1, type1error_interim_2, type1error_interim_3, power_interim_1, power_interim_2, power_interim_3),
               names_to = "effect_sizes", values_to = "value")


newvar_data = str_split(newdata$effect_sizes, pattern = "_", n = nrow(newdata))
newvar_final = str_split(newfinal$effect_sizes, pattern = "_", n = nrow(newfinal))
newvar_interim = str_split(newinterim$effect_sizes, pattern = "_", n = nrow(newinterim))

for(i in 1:nrow(newdata)){
  newdata$stat[i] = newvar_data[[i]][1]
  newdata$stage[i] = newvar_data[[i]][2]
  newdata$effect_sizes[i] = newvar_data[[i]][3]
}

for(i in 1:nrow(newfinal)){
  newfinal$stat[i] = newvar_final[[i]][1]
  newfinal$effect_sizes[i] = newvar_final[[i]][3]
}

for(i in 1:nrow(newinterim)){
  newinterim$stat[i] = newvar_interim[[i]][1]
  newinterim$effect_sizes[i] = newvar_interim[[i]][3]
}

newfinal$effect_sizes = ifelse(newfinal$effect_sizes == 1 & newfinal$stat == "type1error", "0.05/0.05", 
                               ifelse(newfinal$effect_sizes == 1 & newfinal$stat == "power", "0.05/0.20",
                                      ifelse(newfinal$effect_sizes == 2 & newfinal$stat == "type1error", "0.10/0.10", 
                                             ifelse(newfinal$effect_sizes == 2 & newfinal$stat == "power", "0.10/0.25",
                                                    ifelse(newfinal$effect_sizes == 3 & newfinal$stat == "type1error", "0.15/0.15", 
                                                           "0.15/0.30")))))

newinterim$effect_sizes = ifelse(newinterim$effect_sizes == 1 & newinterim$stat == "type1error", "0.05/0.05", 
                               ifelse(newinterim$effect_sizes == 1 & newinterim$stat == "power", "0.05/0.20",
                                      ifelse(newinterim$effect_sizes == 2 & newinterim$stat == "type1error", "0.10/0.10", 
                                             ifelse(newinterim$effect_sizes == 2 & newinterim$stat == "power", "0.10/0.25",
                                                    ifelse(newinterim$effect_sizes == 3 & newinterim$stat == "type1error", "0.15/0.15", 
                                                           "0.15/0.30")))))

####Plots####
library(tidyverse)


#Power vs Sample size for interim analysis
ggplot(data = newinterim[newinterim$stat=="power",], aes(x = n, y = value, color = effect_sizes)) + geom_line(size = 1) + scale_fill_discrete(name = "Effect Sizes") + labs(title = "Power vs Sample Size grouped by Effect Sizes (Interim Analysis)", y = "Power") + theme_bw() + geom_hline(aes(yintercept = 0.8), size = 1)

#Type I error vs Sample size for interim analysis
ggplot(data = newinterim[newinterim$stat=="type1error",], aes(x = n, y = value, color = effect_sizes)) + geom_line(size = 1) + scale_fill_discrete(name = "Effect Sizes") + labs(title = "Type I Error vs Sample Size grouped by Effect Sizes(Interim Analysis)", y = "Type I Error") + theme_bw() + geom_hline(aes(yintercept = 0.005), size = 0.75)


#Power vs Sample size for final analysis
ggplot(data = newfinal[newfinal$stat=="power",], aes(x = n, y = value, color = effect_sizes)) + geom_line(size = 1) + scale_fill_discrete(name = "Effect Sizes") + labs(title = "Power vs Sample Size grouped by Effect Sizes (Final Analysis)", y = "Power") + theme_bw() + geom_hline(aes(yintercept = 0.8), size = 1)

which(newfinal$stat == "power" & newfinal$value >= .80 & newfinal$effect_sizes == "0.15/0.30")

#Type I error vs Sample size for final analysis
ggplot(data = newfinal[newfinal$stat=="type1error",], aes(x = n, y = value, color = effect_sizes)) + geom_line(size = 1) + scale_fill_discrete(name = "Effect Sizes") + labs(title = "Type I Error vs Sample Size grouped by Effect Sizes (Final Analysis)", y = "Type I Error") + theme_bw() + geom_hline(aes(yintercept = 0.048), size = 1)


#Recommended sample size is 122 patients per arm