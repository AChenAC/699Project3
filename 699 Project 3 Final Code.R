library(tidyverse)

###Number of simulations###

(1.96/0.005)^2*(0.2*0.8)

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
  z_stat <-
    (all_data_A - all_data_B) / 
    sqrt(estimated_sigmasq_A / n_A + estimated_sigmasq_B / n_B);
  list(pvalues_t = 2 * pnorm(q = abs(z_stat), mean = 0, sd = 1, 
                          lower.tail = FALSE), 
       runtime = Sys.time() - start_sim);  
}

trial_sim_bernoulli_new <- function(n_sim, 
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
  z_stat <-
    (all_data_A - all_data_B) / 
    sqrt(estimated_sigmasq_A / n_A + estimated_sigmasq_B / n_B);
  pvalues_t = 2 * pnorm(q = abs(z_stat), 
                          mean = 0, sd = 1, 
                          lower.tail = FALSE);
  return(data.frame(z_stat, pvalues_t))
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

for(i in 55:65){
  
  pvals_interim_null_1[[i]] = unname(unlist(trial_sim_bernoulli(25000, i, i, 0.05, 0.05, seed = 456)))
  pvals_interim_null_2[[i]] = unname(unlist(trial_sim_bernoulli(25000, i, i, 0.1, 0.1, seed = 456)))
  pvals_interim_null_3[[i]] = unname(unlist(trial_sim_bernoulli(25000, i, i, 0.15, 0.15, seed = 456)))
  
  type1error_interim_1[i] = length(which(pvals_interim_null_1[[i]]<0.005))/length(pvals_interim_null_1[[i]])
  type1error_interim_2[i] = length(which(pvals_interim_null_2[[i]]<0.005))/length(pvals_interim_null_2[[i]])
  type1error_interim_3[i] = length(which(pvals_interim_null_3[[i]]<0.005))/length(pvals_interim_null_3[[i]])
  
  pvals_interim_alt_1[[i]] = unname(unlist(trial_sim_bernoulli(25000, i, i, 0.05, 0.2, seed = 456)))
  pvals_interim_alt_2[[i]] = unname(unlist(trial_sim_bernoulli(25000, i, i, 0.1, 0.25, seed = 456)))
  pvals_interim_alt_3[[i]] = unname(unlist(trial_sim_bernoulli(25000, i, i, 0.15, 0.3, seed = 456)))
  
  power_interim_1[i] = length(which(pvals_interim_alt_1[[i]]<0.005))/length(pvals_interim_alt_1[[i]])
  power_interim_2[i] = length(which(pvals_interim_alt_2[[i]]<0.005))/length(pvals_interim_alt_2[[i]])
  power_interim_3[i] = length(which(pvals_interim_alt_3[[i]]<0.005))/length(pvals_interim_alt_3[[i]])
}

type1error_interim_1 = type1error_interim_1[55:65]
type1error_interim_2 = type1error_interim_2[55:65]
type1error_interim_3 = type1error_interim_3[55:65]

power_interim_1 = power_interim_1[55:65]
power_interim_2 = power_interim_2[55:65]
power_interim_3 = power_interim_3[55:65]


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
  
  pvals_final_null_1[[i]] = unname(unlist(trial_sim_bernoulli(25000, i, i, 0.05, 0.05, seed = 456)))
  pvals_final_null_2[[i]] = unname(unlist(trial_sim_bernoulli(25000, i, i, 0.1, 0.1, seed = 456)))
  pvals_final_null_3[[i]] = unname(unlist(trial_sim_bernoulli(25000, i, i, 0.15, 0.15, seed = 456)))
  
  type1error_final_1[i] = length(which(pvals_final_null_1[[i]]<0.048))/length(pvals_final_null_1[[i]])
  type1error_final_2[i] = length(which(pvals_final_null_2[[i]]<0.048))/length(pvals_final_null_2[[i]])
  type1error_final_3[i] = length(which(pvals_final_null_3[[i]]<0.048))/length(pvals_final_null_3[[i]])
  
  pvals_final_alt_1[[i]] = unname(unlist(trial_sim_bernoulli(25000, i, i, 0.05, 0.2, seed = 456)))
  pvals_final_alt_2[[i]] = unname(unlist(trial_sim_bernoulli(25000, i, i, 0.1, 0.25, seed = 456)))
  pvals_final_alt_3[[i]] = unname(unlist(trial_sim_bernoulli(25000, i, i, 0.15, 0.3, seed = 456)))
  
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

final = cbind(type1error_final, power_final)
interim = cbind(type1error_interim, power_interim)

final$n = 110:130
interim$n = 55:65

newfinal = final %>%
  pivot_longer(c(type1error_final_1, type1error_final_2, type1error_final_3, power_final_1, power_final_2, power_final_3),
               names_to = "effect_sizes", values_to = "value")

newinterim = interim %>%
  pivot_longer(c(type1error_interim_1, type1error_interim_2, type1error_interim_3, power_interim_1, power_interim_2, power_interim_3),
               names_to = "effect_sizes", values_to = "value")


newvar_final = str_split(newfinal$effect_sizes, pattern = "_", n = nrow(newfinal))
newvar_interim = str_split(newinterim$effect_sizes, pattern = "_", n = nrow(newinterim))

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
ggplot(data = newinterim[newinterim$stat=="power",], aes(x = n, y = value, color = effect_sizes)) + geom_line(size = 1) + scale_fill_discrete(name = "Effect Sizes") + labs(title = "Power vs Sample Size grouped by Effect Sizes (Interim Analysis)", y = "Power") + theme_bw()

#Type I error vs Sample size for interim analysis
ggplot(data = newinterim[newinterim$stat=="type1error",], aes(x = n, y = value, color = effect_sizes)) + geom_line(size = 1) + scale_fill_discrete(name = "Effect Sizes") + labs(title = "Type I Error vs Sample Size grouped by Effect Sizes(Interim Analysis)", y = "Type I Error") + theme_bw() + geom_hline(aes(yintercept = 0.005), size = 0.75)


#Power vs Sample size for final analysis
ggplot(data = newfinal[newfinal$stat=="power",], aes(x = n, y = value, color = effect_sizes)) + geom_line(size = 1) + scale_fill_discrete(name = "Effect Sizes") + labs(title = "Power vs Sample Size grouped by Effect Sizes (Final Analysis)", y = "Power") + theme_bw() + geom_hline(aes(yintercept = 0.8), size = 1)

which(newfinal$stat == "power" & newfinal$value >= .80 & newfinal$effect_sizes == "0.15/0.30")

#Type I error vs Sample size for final analysis
ggplot(data = newfinal[newfinal$stat=="type1error",], aes(x = n, y = value, color = effect_sizes)) + geom_line(size = 1) + scale_fill_discrete(name = "Effect Sizes") + labs(title = "Type I Error vs Sample Size grouped by Effect Sizes (Final Analysis)", y = "Type I Error") + theme_bw() + geom_hline(aes(yintercept = 0.048), size = 1)



#Facet Wrapping

newfinal$stage = "final"
newinterim$stage = "interim"

newdata = rbind(newfinal, newinterim)

ggplot(data = newdata, aes(x = n, y = value, color = effect_sizes)) + geom_line(size = 1.1) + theme_bw() + facet_grid(stat ~ stage, scales = "free", ) + theme(text = element_text(size = 20))

#Recommended sample size is 120 patients per arm


###Probability of each outcome###

df_interim_null_1 = trial_sim_bernoulli_new(25000, 60, 60, 0.05, 0.05, seed = 456)
df_interim_null_2 = trial_sim_bernoulli_new(25000, 60, 60, 0.1, 0.1, seed = 456)
df_interim_null_3 = trial_sim_bernoulli_new(25000, 60, 60, 0.15, 0.15, seed = 456)
  
zstat_interim_null_1 = df_interim_null_1$z_stat
zstat_interim_null_2 = df_interim_null_2$z_stat
zstat_interim_null_3 = df_interim_null_3$z_stat
  
df_interim_alt_1 = trial_sim_bernoulli_new(25000, 60, 60, 0.05, 0.2, seed = 456)
df_interim_alt_2 = trial_sim_bernoulli_new(25000, 60, 60, 0.1, 0.25, seed = 456)
df_interim_alt_3 = trial_sim_bernoulli_new(25000, 60, 60, 0.15, 0.3, seed = 456)

zstat_interim_alt_1 = df_interim_alt_1$z_stat
zstat_interim_alt_2 = df_interim_alt_2$z_stat
zstat_interim_alt_3 = df_interim_alt_3$z_stat


#These are the probabilities of stopping at the interim under the null
prop_SOC_better_interim_null_1 = sum(zstat_interim_null_1>2.81, na.rm = T)/length(zstat_interim_null_1) #Removing 60 missing values
prop_trt_better_interim_null_1 = sum(zstat_interim_null_1 < -2.81, na.rm = T)/length(zstat_interim_null_1) #Removing 60 missing values

prop_SOC_better_interim_null_2 = sum(zstat_interim_null_2>2.81)/length(zstat_interim_null_2)
prop_trt_better_interim_null_2 = sum(zstat_interim_null_2 < -2.81)/length(zstat_interim_null_2)

prop_SOC_better_interim_null_3 = sum(zstat_interim_null_3>2.81)/length(zstat_interim_null_3)
prop_trt_better_interim_null_3 = sum(zstat_interim_null_3 < -2.81)/length(zstat_interim_null_3)

#These are the probabilities of stopping at the interim under the alternative
prop_SOC_better_interim_alt_1 = sum(zstat_interim_alt_1>2.81)/length(zstat_interim_alt_1)
prop_trt_better_interim_alt_1 = sum(zstat_interim_alt_1 < -2.81)/length(zstat_interim_alt_1)

prop_SOC_better_interim_alt_2 = sum(zstat_interim_alt_2>2.81)/length(zstat_interim_alt_2)
prop_trt_better_interim_alt_2 = sum(zstat_interim_alt_2 < -2.81)/length(zstat_interim_alt_2)

prop_SOC_better_interim_alt_3 = sum(zstat_interim_alt_3>2.81)/length(zstat_interim_alt_3)
prop_trt_better_interim_alt_3 = sum(zstat_interim_alt_3 < -2.81)/length(zstat_interim_alt_3)


#Final Analysis

df_final_null_1 = trial_sim_bernoulli_new(25000, 120, 120, 0.05, 0.05, seed = 456)
df_final_null_2 = trial_sim_bernoulli_new(25000, 120, 120, 0.1, 0.1, seed = 456)
df_final_null_3 = trial_sim_bernoulli_new(25000, 120, 120, 0.15, 0.15, seed = 456)

zstat_final_null_1 = df_final_null_1$z_stat
zstat_final_null_2 = df_final_null_2$z_stat
zstat_final_null_3 = df_final_null_3$z_stat

df_final_alt_1 = trial_sim_bernoulli_new(25000, 120, 120, 0.05, 0.2, seed = 456)
df_final_alt_2 = trial_sim_bernoulli_new(25000, 120, 120, 0.1, 0.25, seed = 456)
df_final_alt_3 = trial_sim_bernoulli_new(25000, 120, 120, 0.15, 0.3, seed = 456)

zstat_final_alt_1 = df_final_alt_1$z_stat
zstat_final_alt_2 = df_final_alt_2$z_stat
zstat_final_alt_3 = df_final_alt_3$z_stat


#These are the probabilities of rejecting the null (there is a treatment difference) at the end of the trial 
#under the null
prop_SOC_better_final_null_1 = sum(zstat_final_null_1>1.98)/length(zstat_final_null_1)
prop_trt_better_final_null_1 = sum(zstat_final_null_1 < -1.98)/length(zstat_final_null_1)

prop_SOC_better_final_null_2 = sum(zstat_final_null_2>1.98)/length(zstat_final_null_2)
prop_trt_better_final_null_2 = sum(zstat_final_null_2 < -1.98)/length(zstat_final_null_2)

prop_SOC_better_final_null_3 = sum(zstat_final_null_3>1.98)/length(zstat_final_null_3)
prop_trt_better_final_null_3 = sum(zstat_final_null_3 < -1.98)/length(zstat_final_null_3)


#These are the probabilities of rejecting the null (there is a treatment difference) at the end of the trial
#under the alternative.
prop_SOC_better_final_alt_1 = sum(zstat_final_alt_1>1.98)/length(zstat_final_alt_1)
prop_trt_better_final_alt_1 = sum(zstat_final_alt_1 < -1.98)/length(zstat_final_alt_1)

prop_SOC_better_final_alt_2 = sum(zstat_final_alt_2>1.98)/length(zstat_final_alt_2)
prop_trt_better_final_alt_2 = sum(zstat_final_alt_2 < -1.98)/length(zstat_final_alt_2)

prop_SOC_better_final_alt_3 = sum(zstat_final_alt_3>1.98)/length(zstat_final_alt_3)
prop_trt_better_final_alt_3 = sum(zstat_final_alt_3 < -1.98)/length(zstat_final_alt_3)


#This is the probability of failing to reject the null (i.e. no difference) at the end of the trial under the null

prop_fail_to_reject_null_1 = 1 - prop_SOC_better_final_null_1 - prop_trt_better_final_null_1
prop_fail_to_reject_null_2 = 1 - prop_SOC_better_final_null_2 - prop_trt_better_final_null_2
prop_fail_to_reject_null_3 = 1 - prop_SOC_better_final_null_3 - prop_trt_better_final_null_3

#This is the probability of failing to reject the null (i.e. no difference) at the end of the trial under the alternative

prop_fail_to_reject_alt_1 = 1 - prop_SOC_better_final_alt_1 - prop_trt_better_final_alt_1
prop_fail_to_reject_alt_2 = 1 - prop_SOC_better_final_alt_2 - prop_trt_better_final_alt_2
prop_fail_to_reject_alt_3 = 1 - prop_SOC_better_final_alt_3 - prop_trt_better_final_alt_3



###Expected Sample Size###

expected_1 = 60*prop_trt_better_interim_alt_1 + 120*(1-prop_trt_better_interim_alt_1)
expected_2 = 60*prop_trt_better_interim_alt_2 + 120*(1-prop_trt_better_interim_alt_2)
expected_3 = 60*prop_trt_better_interim_alt_3 + 120*(1-prop_trt_better_interim_alt_3)

se_expected_1 = 60*sqrt((prop_trt_better_interim_alt_1*(1-prop_trt_better_interim_alt_1))/25000)
se_expected_2 = 60*sqrt((prop_trt_better_interim_alt_2*(1-prop_trt_better_interim_alt_2))/25000)
se_expected_3 = 60*sqrt((prop_trt_better_interim_alt_3*(1-prop_trt_better_interim_alt_3))/25000)

var_expected_1 = se_expected_1^2
var_expected_2 = se_expected_2^2
var_expected_3 = se_expected_3^2