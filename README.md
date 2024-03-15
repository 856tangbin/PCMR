# PCMR: Pleiotropic Clustering of Mendelian Randomization
## Installation

```
# install.packages("devtools")
library(devtools)
install_github("856tangbin/PCMR")
```





```
# install.packages("C:/Users/BinTang/Desktop/PCMR/PCMR_0.1.0.tar.gz", repos = NULL, type = "source")
library(PCMR)

# read data,true gamma:0 and 0.707. And proportion are 1:1.
X_clump = read.csv("IVs.csv")
X_clump1 = read.csv("initEst.csv")

# Part 1. initial value of sigma2 and rho
init = PCMR_initEst(X_clump1$beta_hat_1,X_clump1$seb1,
             X_clump1$beta_hat_2,X_clump1$seb2)

init

# Part 2. fitting
# The fixed effect model, model = "2" , simplified
result_fixed = PCMR(X_clump$beta_hat_1, X_clump$seb1,
              X_clump$beta_hat_2,X_clump$seb2,num_gamma = 2,model="2",
              sigma2 = init$sigma2)

# Causal testing
result_fixed = PCMR_testCausal(result_fixed)

print(c(result_fixed$gamma,result_fixed$pi_gamma,result_fixed$Pvalue))
# Plotting
PCMR_plot(result_fixed)


# The random effect model, model = "2" , intact
result_random = PCMR(X_clump$beta_hat_1, X_clump$seb1,
              X_clump$beta_hat_2,X_clump$seb2,num_gamma = 2,model="1",
              isIntact=T,rho=init$rho,sigma2 = init$sigma2)

# Causal testing
result_random = PCMR_testCausal(result_random)

print(c(result_random$gamma,result_random$pi_gamma,result_random$Pvalue))
# Plotting
PCMR_plot(result_random)


# Part 3*(可选). Set initial parameters as required 

# setting the max gSigma2 ####################################################
result0 = PCMR(X_clump$beta_hat_1, X_clump$seb1,
               X_clump$beta_hat_2,X_clump$seb2,num_gamma = 1,model="1",isIntact=T,rho=init$rho,
               sigma2 = init$sigma2)

result1 = PCMR(X_clump$beta_hat_1, X_clump$seb1,
                     X_clump$beta_hat_2,X_clump$seb2,num_gamma = 2,model="1",
                     isIntact=T,rho=init$rho,sigma2 = init$sigma2,
                     max_gSigma2 = result0$gSigma2)

PCMR_plot(result1)
# setting initial values 0 and 0.01 for gamma.
# Note: Different initial values may lead to local convergence without affecting the final convergence result.
result1 = PCMR(X_clump$beta_hat_1, X_clump$seb1,
               X_clump$beta_hat_2,X_clump$seb2,num_gamma = 2,model="1",
               isIntact=T,rho=init$rho,sigma2 = init$sigma2,
               gamma=c(0,0.01))
PCMR_plot(result1)

# setting initial values 0.1 for b.
result1 = PCMR(X_clump$beta_hat_1, X_clump$seb1,
               X_clump$beta_hat_2,X_clump$seb2,num_gamma = 2,model="1",
               isIntact=T,rho=init$rho,sigma2 = init$sigma2,
               b=0.1)
PCMR_plot(result1)

# Part 4. Detecting the presence of correlated horizontal pleiotropy
library(parallel)
result_random = PCMR_testCorPlei(result_random,samples = 100)

print(result_random$CHVP_test)

# correct P-value of PCMR's pleiotropy test
result_random = PCMR_correct(result_random,ref_beta_outcome = X_clump1$beta_hat_2,ref_se_outcome = X_clump1$seb2,samples = 20, 
             sample_boot = 30)

print(result_random$CHVP_test_correct) # corrected P-value in detecting correlated horizontal pleiotropy.
print(result_random$CHVP_test_correctRange) # Error range of the P-value 




```

