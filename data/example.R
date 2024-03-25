library(PCMR)
library(parallel)

X_clump = read.table("./data/scz_mdd/IVs_scz_mdd.csv",
                     header = 1,
                     sep = ",")
X_clump1 = read.table("./data/scz_mdd/initEst_scz_mdd.csv",
                      header = 1,
                      sep = ",")
set.seed(100)

# Part 1 estimate initial values
init = PCMR_initEst(X_clump1$beta_hat_1,X_clump1$seb1,
                    X_clump1$beta_hat_2,X_clump1$seb2)

# Part 2 model fit
result_random = PCMR(X_clump$beta_hat_1, X_clump$seb1,
                     X_clump$beta_hat_2,X_clump$seb2,num_gamma = 2,model="1",
                     isIntact=T,rho=init$rho,sigma2 = init$sigma2) # random effect model
result_random$gamma
result_random$pi_gamma
PCMR_plot(result_random)

# Part 2 Heterogeneity test in detecting correlated horizontal pleiotropy
result_random = PCMR_cEst(result_random,
                         ref_beta_outcome = X_clump1$beta_hat_2,ref_se_outcome = X_clump1$seb2,
                         samples=20,sample_boot = 30,cores = 8) # estimate the factor c
print(c(result_random$c,result_random$c_sd))

result_random = PCMR_testCausal_bootstrap(result_random,samples=100,cores=8) # bootstrapping to estimate D_HVP
print(c(result_random$D_HVP))

result_random = PCMR_testCorPlei(result_random)
print(c(result_random$CHVP_test,result_random$CHVP_test_Range))


# Part 3 Causality analysis

result_random = PCMR_testCausal(result_random)
print(result_random$Pvalue)


# Part 4 classifed IV categories

prb_thrd = 0.5
probability_1 = rowSums(result_random$W[,1,]) # the probability of IVs belonging to 1th IV category

stringi::stri_c(X_clump$rsid[probability_1 > prb_thrd],collapse = " ")
stringi::stri_c(X_clump$rsid[probability_1 < 1 - prb_thrd],collapse = " ")


