library(TwoSampleMR)

library(PCMR)
library(parallel)

opengwas_jwt = "" # Token (JWT) by registered on https://api.opengwas.io/profile/

set.seed(0)
# background variants information ################
# We recommend using variants with P-values greater than 0.5 in the GWAS of exposure and performing LD prunning as a background for assessing uncorrelated pleiotropy.
data(IVs_and_InitEst,package="PCMR") # Due to network limitations, we randomly screen a small number of genetic variants for background display here.
exposure_dat1 <- extract_outcome_data(snps=sample(X_clump1$rsid,500), outcomes="ieu-b-4877",opengwas_jwt = opengwas_jwt) # smoking
exposure_dat1 = exposure_dat1[exposure_dat1$pval.outcome > 0.5,]
outcome_dat1 <- extract_outcome_data(snps=exposure_dat1$SNP, outcomes="ieu-b-5110",opengwas_jwt = opengwas_jwt) # BIP
dat1 = merge(exposure_dat1,outcome_dat1,by = "SNP")


X_clump1 = data.frame(rsid = dat1$SNP,
                      beta_hat_1 = dat1$beta.outcome.x,
                      seb1 = dat1$se.outcome.x,
                      beta_hat_2 = dat1$beta.outcome.y,
                      seb2 = dat1$se.outcome.y)

# IVs information ##################
exposure_dat <- extract_instruments("ieu-b-4877",r2 = 0.1,opengwas_jwt = opengwas_jwt) # smoking
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-b-5110",opengwas_jwt = opengwas_jwt) # BIP

dat <- harmonise_data(exposure_dat, outcome_dat)
X_clump = data.frame(rsid = dat$SNP,
                     beta_hat_1 = dat$beta.exposure,
                     seb1 = dat$se.exposure,
                     beta_hat_2 = dat$beta.outcome,
                     seb2 = dat$se.outcome)

# Perform PCMR ###################
# 0. Estimate initial values.
init = PCMR_initEst(X_clump1$beta_hat_1,X_clump1$seb1,
                    X_clump1$beta_hat_2,X_clump1$seb2)

# 1. Pleiotropic clustering;
result_random = PCMR(X_clump$beta_hat_1, X_clump$seb1,
                     X_clump$beta_hat_2,X_clump$seb2,num_gamma = 2,model="1",
                     isIntact=T,rho=init$rho,sigma2 = init$sigma2)

PCMR_plot(result_random)

# 2. Heterogeneity test in detecting correlated horizontal pleiotropy.
result_random = PCMR_cEst(result_random,ref_beta_outcome = X_clump1$beta_hat_2,ref_se_outcome = X_clump1$seb2,cores=10) # Estimate the factor c.
result_random = PCMR_testCausal_bootstrap(result_random,cores=10) # Bootstrapping to estimate D_HVP.
result_random = PCMR_testCorPlei(result_random) # Calculate Pvalue of heterogeneity test according to c and D_HVP.

# 3. Causality evaluation in the presence of correlated horizontal pleiotropy.
# Recommended to be applied in the presence of correlated horizontal pleiotropy, e.g. P_{plei-test} <= 0.20.
result_random = PCMR_testCausal(result_random)
