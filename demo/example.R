library(PCMR)
library(parallel)

t = Sys.time()
X_clump = read.table("./data/scz_mdd/IVs_scz_mdd.csv",
                     header = 1,
                     sep = ",")
X_clump1 = read.table("./data/scz_mdd/initEst_scz_mdd.csv",
                      header = 1,
                      sep = ",")

set.seed(0)

# Part 0 estimate initial values
init = PCMR_initEst(X_clump1$beta_hat_1,X_clump1$seb1,
                    X_clump1$beta_hat_2,X_clump1$seb2)

# Part 1 model fit
result_random = PCMR(X_clump$beta_hat_1, X_clump$seb1,
                     X_clump$beta_hat_2,X_clump$seb2,num_gamma = 2,model="1",
                     isIntact=T,rho=init$rho,sigma2 = init$sigma2) # random effect model
result_random$gamma
result_random$pi_gamma

# Part 2 Heterogeneity test in detecting correlated horizontal pleiotropy
result_random = PCMR_cEst(result_random,ref_beta_outcome = X_clump1$beta_hat_2,ref_se_outcome = X_clump1$seb2,cores=10) # estimate the factor c
result_random = PCMR_testCausal_bootstrap(result_random,cores=10) # bootstrapping to estimate D_HVP
result_random = PCMR_testCorPlei(result_random) # calculate Pvalue of heterogeneity


# Part 3 Causality analysis
result_random = PCMR_testCausal(result_random)
print(Sys.time() - t)


# Part 4 classifed IV categories

prb_thrd = 0.5
probability_1 = rowSums(result_random$W[,1,]) # the probability of IVs belonging to 1th IV category

stringi::stri_c(X_clump$rsid[probability_1 > prb_thrd],collapse = " ")
stringi::stri_c(X_clump$rsid[probability_1 < 1 - prb_thrd],collapse = " ")


