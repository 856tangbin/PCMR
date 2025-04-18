# PCMR: Pleiotropic Clustering of Mendelian Randomization

> This R package implements the PCMR method (https://doi.org/10.1038/s41467-025-57912-5) and provides efficient solutions and new insights for MR in clinically relevant or genetically overlapping traits. 

## Usage 

PCMR is a clustering model for addressing various horizontal or vertical pleiotropic (HVP) effects. PCMR contains three components: 

- **Pleiotropic Clustering** (PCMR): Cluster instruments using the function of `PCMR`;
- **Detecting correlated horizontal pleiotropy** (PCMR's pleiotropy test): Detect the presence of correlated horizontal pleiotropy using the function of `PCMR_testCorPlei`;
- **Causal analysis** (PCMR's causality evaluation): Evaluate whether a discernible dominant IV category supports a non-zero causal effect using the function of `PCMR_testCausal`. Recommended to be applied in the presence of correlated horizontal pleiotropy, e.g. $P_{plei-test} <= 0.20$. 

Data preprocessing of IVs for PCMR was provided in `./demo/IVs_filter.R`; We provided an example of SCZ and MDD in `./demo/example.R`, as well as, an example of combining R package `TwoSampleMR` to analysis the causal relationship between Smoking and BIP in `./demo/example_TwoSampleMR.R`.



## Installation

```R
# install.packages("devtools")
library(devtools)
install_github("856tangbin/PCMR")
```





## Example

> In the application from schizophrenia (SCZ) to major depressive disorder (MDD)

Instrument variables: [IVs_scz_mdd.csv](data\scz_mdd\IVs_scz_mdd.csv); Random sample variants: [initEst_scz_mdd.csv](data\scz_mdd\initEst_scz_mdd.csv). 

These instruments and random sample variants are obtained through functions in the R package cause, see`./demo/IVs_filter.R`. 

```R
library(PCMR)
library(parallel)

# X_clump = read.table("./data/scz_mdd/IVs_scz_mdd.csv",header=1,sep=",")
# X_clump1 = read.table("./data/scz_mdd/initEst_scz_mdd.csv",header=1,sep=",")
data(IVs_and_InitEst,package="PCMR")

set.seed(0)
# 0. Estimate initial values.
init = PCMR_initEst(X_clump1$beta_hat_1,X_clump1$seb1,
                    X_clump1$beta_hat_2,X_clump1$seb2)

# 1. Pleiotropic clustering; 
result_random = PCMR(X_clump$beta_hat_1, X_clump$seb1,
              X_clump$beta_hat_2,X_clump$seb2,num_gamma = 2,model="1",
              isIntact=T,rho=init$rho,sigma2 = init$sigma2)

# 2. Heterogeneity test in detecting correlated horizontal pleiotropy.
result_random = PCMR_cEst(result_random,ref_beta_outcome = X_clump1$beta_hat_2,ref_se_outcome = X_clump1$seb2,cores=10) # Estimate the factor c.
result_random = PCMR_testCausal_bootstrap(result_random,cores=10) # Bootstrapping to estimate D_HVP.
result_random = PCMR_testCorPlei(result_random) # Calculate Pvalue of heterogeneity test according to c and D_HVP.

# 3. Causality evaluation in the presence of correlated horizontal pleiotropy.
# Recommended to be applied in the presence of correlated horizontal pleiotropy, e.g. P_{plei-test} <= 0.20. 
result_random = PCMR_testCausal(result_random)
```

*Note: The Step 0, Step 1 and Step 3 were quick, about one minute in total. The Step 2 for detecting correlated horizontal pleiotropy took about 30+ minutes. You can also directly load the results of the completed calculations using* `data(scz_mdd_results,package="PCMR")`.





## Results

### 0. Estimating initial value

```
> init
$sigma2
[1] 0.000000e+00 1.992793e-05

$rho
[1] 0.01972954
```

***Note***: The estimates of `sigma2` and `rho` represent the initial values of uncorrelated horizontal pleiotropy and sample overlap, respectively. 



### 1. Pleiotropic Clustering

```R
> result_random$gamma # Correlated HVP effects
 33.33333%  66.66667% 
0.01703264 0.18722993 

> result_random$gSigma2
[1] 0.004008198 0.002144189

> result_random$pi_gamma
[1] 0.4981254 0.5018746

> PCMR_plot(result_random)
```

***Note***: The classified estimates of `gamma` are the sum of correlated horizontal pleiotropic effect and causal effect, being Correlated HVP effect; `gSigma2` and `pi_gamma` measure the random correlated HVP effects and the proportion of distinct categories, respectively. The default model is the random effect model of PCMR (`model=1`), while the fixed model of PCMR (`model=2`) sets `gSigma2` at zero. 

<img src="README.assets/Rplot_random.svg" alt="Rplot_random" align=center />



### 2. Heterogeneity test in detecting correlated horizontal pleiotropy

```R
> print(result_random$CHVP_test) # Pvalue of heterogeniety test
[1] 1.128664e-07
```

***Note***: The heterogeneity test implies that there is insignificantly correlated horizontal pleiotropy between SCZ and MDD ($P_{plei-test}=9.00\times 10^{-08}$). 



### 3.  Causality evaluation in the presence of correlated horizontal pleiotropy

```R
> print(result_random$Pvalue)
[1] 0.3290909

> print(c(result_random$effect,result_random$discernable_prob))
0.1872299 0.7917918 
```

***Note***: `Pvalue` is the statistic value of PCMR's causality evaluation, which evaluates relationships in the presence of correlated horizontal pleiotropy ($P_{plei-test} < 0.20$).  The value of `effect` indicates the estimation of causal effect supported by the largest group (BLUE category), and `discernable_prob` represents the probability that the largest sample IV group is discernable as the dominant population IV group. 



### 4. Conclusions

In SCZ and MDD, 

- **Pleiotropic Clustering** (PCMR): PCMR clusters the instrumental variables into two distinct categories, with correlated HVP effects of 0.017 and 0.187;
- **Heterogeneity test** (PCMR's pleiotropy test): The difference of correlated HVP effects is significant ($P_{plei-test}=1.13\times 10^{-07}$), meaning the presence of correlated horizontal pleiotropy. 
- **Causal analysis** (PCMR's causality evaluation): The results of PCMR's causality evaluation indicated that the relationship of SCZ on MDD was insignificant ($P=0.329$). The estimated causal effect supported by the largest IV group is 0.187, but the discernable probability is solely 0.791. 





## Extension: Integrating biological information for enhancing causal inference

The classified IV categories by PCMR facilitate the integration of biological information for mechanism interpretation, offering a valuable avenue to exclude correlated horizontal pleiotropic variants for enhancing causal inference. For example: 

The instrument classified by PCMR can be mapped into genes at the website: https://biit.cs.ut.ee/gprofiler/snpense, and then using those genes for enrichment analysis of biological process at the website: https://biit.cs.ut.ee/gprofiler/gost. The enrichment analysis may aid in excluding correlated horizontal pleiotropic variants to enhance causal inference.

```R
prb_thrd = 0.5
probability_1 = rowSums(result_random$W[,1,]) # the probability of IVs belonging to 1th IV category

stringi::stri_c(X_clump$rsid[probability_1 > prb_thrd],collapse = " ") # The IV category with correlated HVP effect of result$gamma[1] 
stringi::stri_c(X_clump$rsid[probability_1 < 1 - prb_thrd],collapse = " ") # The IV category with correlated HVP effect of result$gamma[2] 
```

**The IV category with correlated HVP effect of 0.01703264 (GRAY)**

```R
> stringi::stri_c(X_clump$rsid[probability_1 > prb_thrd],collapse = " ")
[1] "rs113113059 rs8073146 rs2840676 rs7515363 rs56335113 rs1915019 rs308697 rs1615350 rs6001259 rs167924 rs9876421 rs6549963 rs6538539 rs7575796 rs10086619 rs11210892 rs77492327 rs217336 rs2381411 rs39967 rs2910032 rs72802868 rs79212538 rs12771371 rs6984242 rs58120505 rs17731 rs4766428 rs16940992 rs10415576 rs2999392 rs3791710 rs2252074 rs12327967 rs778371 rs4575535 rs10117 rs9687282 rs416571 rs9393698 rs764284 rs12540417 rs7801375 rs2894222 rs9270836 rs215412 rs7647398 rs1430894 rs80094991 rs5751191 rs1451488 rs7563610 rs56245805 rs12991836 rs16825349 rs500102 rs7090337 rs9454727 rs324017 rs61937595 rs2013949 rs1901512 rs10861176 rs2455415 rs10035564 rs34611983 rs17194490 rs12285419 rs3729986 rs6926151 rs71301816 rs4779050 rs867810 rs3795310 rs246024 rs11993663 rs4921741 rs72974238 rs4902961 rs2965189 rs72986630 rs6482437 rs16867571 rs7803571 rs2944839 rs11954859 rs73229090 rs113264400 rs1792709 rs17571951 rs12883788 rs2053079 rs12865628 rs6546857 rs2802535 rs62183855 rs34879738 rs2935244 rs10108980 rs72980087 rs7115692 rs893949 rs7112912 rs3016386 rs79445414 rs713692 rs6798742 rs9813516 rs2909457 rs35734242"
```

![gProfiler_hsapiens_2024-04-24_12-16-52](README.assets/gProfiler_hsapiens_2024-04-24_12-16-52.png)

**The IV category with correlated HVP effect of 0.18722993 (BLUE)**


```R
> stringi::stri_c(X_clump$rsid[probability_1 < 1 - prb_thrd],collapse = " ")
[1] "rs6715366 rs12432904 rs12892189 rs11687313 rs3739118 rs11263861 rs12310367 rs10777956 rs4702 rs67439964 rs11136325 rs4129585 rs3824451 rs7927176 rs12293670 rs12652777 rs149165 rs7830315 rs11941714 rs13145415 rs117145318 rs2333321 rs61405217 rs12498839 rs6943762 rs13233308 rs16851048 rs4636654 rs11027839 rs2456020 rs198806 rs12199613 rs1796518 rs6909479 rs4712938 rs2531806 rs13195636 rs35531336 rs1265099 rs9272446 rs9461856 rs116416291 rs9461916 rs570263 rs9274623 rs209474 rs1144708 rs2967 rs11693094 rs12129573 rs6762456 rs12489270 rs5995756 rs9623320 rs699318 rs4812325 rs60135207 rs13016542 rs7900775 rs11191580 rs79780963 rs2815731 rs9910403 rs8055219 rs187557 rs10860960 rs1604060 rs2241033 rs56393513 rs4936277 rs4987094 rs2514218 rs1881046 rs79210963 rs634940 rs9304548 rs2710323 rs13080668 rs7634476 rs11638554 rs708228 rs35351411 rs3814883 rs2190864 rs7191183 rs72692857 rs13107325 rs7403630 rs11807834 rs72769124 rs145071536 rs9428966 rs4658559 rs2489213 rs7099380 rs6974218 rs74480281 rs11682175 rs12969453 rs74914300 rs72934602 rs4632195 rs9636107 rs505061 rs9318627 rs13407231 rs6701322 rs56205728 rs2255663 rs117472063 rs11664298 rs10894308 rs4936215 rs7001340 rs61786047 rs6520064 rs2238057 rs740417 rs12712510 rs6125656 rs926288"
```
![gProfiler_hsapiens_2024-04-24_12-15-08](README.assets/gProfiler_hsapiens_2024-04-24_12-15-08.png)

The BLUE category showed significant enrichment in biological processes primarily related to signaling transmission, with chemical synaptic transmission being the second significant ($2.344\times 10^{-3}$), while the GRAY category exhibited enrichment in two biological processes. As psychiatric disorders are associated with signal transmission, the BLUE category with a larger correlated HVP effect might exhibit correlated horizontal pleiotropy. 

PCMR also contains a test based on bootstrapping against a particular correlated HVP effect. (**If the effect is determined to be the causal effect, the bootstrapping test is a causal inference**)：

```R
# estimate and standard error: 
> print(c(result_random$bootstrap$mean_minClass,result_random$bootstrap$sd_minClass)) # Bootstrapping test for the smaller correlated HVP effect; When replacing `min` with `max` is the test for the larger correlated HVP effect.
[1] 0.01915753 0.02209725

> result_random$bootstrap$P_min
[1] 0.386865
```

**Conclusion**: Based on the smaller correlated HVP effect, the estimates (standard error) is 0.019 (0.022) and SCZ has an insignificant relationship to MDD (P=0.387). 





## More Information



The software has been tested on Windows 10 with 13th Gen Intel(R) Core(TM) i7-13700K  3.40 GHz and R version 4.1.2.
