#' @title Pleiotropy Clustering Mendelian Randomiztion(PCMR)
#' @author bintang<tangbin0513@gmail.com>
#'
#' @param beta_ex A numeric vector of beta-coefficient values for genetic associations with the exposure.
#' @param beta_ex_se The standard errors associated with the beta-coefficients beta_ex.
#' @param beta_out A numeric vector of beta-coefficient values for genetic associations with the outcome.
#' @param beta_out_se The standard errors associated with the beta-coefficients beta_out.
#' @param num_gamma The number of correlated HVP effect, being the sum of correlated horizontal pleiotropic and causal effects. 1 represents assuming IVs only have causal effect. 2 represents assuming IVs exist correlated pleiotropy besides causal effect.
#' @param num_b The number of classes of uncorrelated pleiotropy.
#' @param model A character representing different model. "1" represents PCMR model with random sigma2 and gSigma2. "2" represents PCMR model with random sigma2 and fixed gSigma2."3" represents PCMR model with fixed sigma2 and random gSigma2."4" represents PCMR model with fixed sigma2 and gSigma2. Default is "2".
#' @param isIntact Whether using the intact model to fitting. Default is FALSE.
#' @param rho The covariance of beta_ex and beta_out. Default is zero.
#' @param gamma The initial value of correlated pleiotropy and causal effect.  A numeric vector of same length as num_gamma. The default is quantile.
#' @param gSigma2 The intitial variance of gamma. The default is zero.
#' @param pi_gamma The initial proportion of gamma. The default is the same ratio for each gamma.
#' @param b A numeric value of uncorrelated pleiotropy. The default is zero.
#' @param sigma2 A bivector and the initial variance of b. The default is bivector of zero.
#' @param pi_b The initial proportion of b.The default is the same ratio for each sigma2.
#' @param paraUpdate A list containing which variables can be updated.
#' @param thred The thredhold of convergence contidion. The default is 1e-6.
#' @param maxIteration The max iteration. The default is 1000.
#'
#' @return gamma The estimated of correlated pleiotropy and causal effect.
#' @return gSigma2 The estimated variance of gamma.
#' @return pi_gamma The estimated proportion of gamma.
#' @return b The estimated of uncorrelated pleiotropy.
#' @return sigma2 The estimated variance of b.
#' @return pi_b The estimated proportion of sigma2.

#'
#' @export
#'
#' @references \url{urlrul}
#'

# testing causal effect
PCMR = function(beta_ex, beta_ex_se, beta_out, beta_out_se,num_gamma,model="2",isIntact=FALSE,
                num_b=2, rho=0,
                gamma=NA, pi_gamma=NA, b=NA, pi_b=NA,sigma2=NA,max_sigma2 = NA,gSigma2 = NA,max_gSigma2=NA,
                paraUpdate=NA,
                thred=1e-6, maxIteration=1000){

    if((rho != 0)&(!isIntact)){
        Warning("The covriance of beta_ex and beta_out is not zero. Please using intact model fitting.")
    }

    Paras = list(beta_ex=beta_ex, beta_ex_se=beta_ex_se, beta_out=beta_out, beta_out_se=beta_out_se,num_gamma=num_gamma,model=model,isIntact=isIntact,
                 num_b=num_b, rho=rho,
                 gamma=gamma, pi_gamma=pi_gamma, b=b, pi_b=pi_b,sigma2=sigma2,max_sigma2 = max_sigma2,gSigma2 = gSigma2,max_gSigma2=max_gSigma2,
                 paraUpdate=paraUpdate,
                 thred=thred, maxIteration=maxIteration)

    if(any(is.na(Paras$paraUpdate))){
        if(model == "1"){ # model 1 : gSigma2 T ; sigma2 T
            Paras$paraUpdate=list(gamma=rep(T,num_gamma),
                                  pi_gamma = rep(T,num_gamma),
                                  b=rep(T,num_b),
                                  pi_b = rep(TRUE,num_b),gSigma2=rep(T,num_gamma),
                                  sigma2 = c(F,T))

        }else if(model == "2"){ # model 2 : gSigma2 F ; sigma2 T
            Paras$paraUpdate=list(gamma=rep(T,num_gamma),
                                  pi_gamma = rep(T,num_gamma),
                                  b=rep(T,num_b),
                                  pi_b = rep(TRUE,num_b),gSigma2=rep(F,num_gamma),
                                  sigma2 = c(F,T))
        }else if(model == "3"){ # model 3 : gSigma2 T ; sigma2 F
            Paras$paraUpdate=list(gamma=rep(T,num_gamma),
                                  pi_gamma = rep(T,num_gamma),
                                  b=rep(T,num_b),
                                  pi_b = rep(TRUE,num_b),gSigma2=rep(T,num_gamma),
                                  sigma2 = c(F,F))
        }else if(model == "4"){ # model 4 : gSigma2 F ; sigma2 F
            Paras$paraUpdate=list(gamma=rep(T,num_gamma),
                                  pi_gamma = rep(T,num_gamma),
                                  b=rep(T,num_b),
                                  pi_b = rep(TRUE,num_b),gSigma2=rep(F,num_gamma),
                                  sigma2 = c(F,F))
        }else{
            # 报错
            stop("Please input parameters correctly.")
        }
    }

    # The clustering of PCMR
    result = pleiClassify(Paras)

    if(length(result$gamma) > 1){
        if(abs(result$gamma[1] - result$gamma[2])) < 1e-3){
            print("Since the estimated correlated HVP effects are identical, there may in absence of correlated pleiotropy.")
        }
    }

    return(result)
}

PCMR_initEst = function(beta_ex, beta_ex_se, beta_out, beta_out_se){

    res = PCMR(beta_ex=beta_ex, beta_ex_se=beta_ex_se, beta_out=beta_out, beta_out_se=beta_out_se,isIntact=FALSE,
                 num_gamma = 1, num_b = 2, gamma = 0 , b = c(0), sigma2 = c(0, 0),
                 paraUpdate = list(
                     gamma = rep(F, 1),
                     pi_gamma = rep(T, 1),
                     b = rep(F, 1),
                     pi_b = rep(T, 2),
                     gSigma2 = F,
                     sigma2 = c(F,TRUE)))



    rho = cov(beta_ex/beta_ex_se,beta_out/beta_out_se)
    return(list(sigma2=res$sigma2,rho=rho))
}


PCMR_testCausal = function(result){
    # save estimated value
    Paras0 = result$Paras
    Paras0$gamma=result$gamma
    Paras0$pi_gamma= result$pi_gamma
    Paras0$b = result$b
    Paras0$pi_b = result$pi_b
    Paras0$sigma2 = result$sigma2
    Paras0$gSigma2 = result$gSigma2

    Paras0$paraUpdate$pi_gamma = rep(FALSE,Paras0$num_gamma) # Assuming the proportion of gamma is true from first estimattion.

    # Selecting the higher proportion of gamma be causal effect. And test H0 by likelihood-ratio test: causal effect whether is zero.
    argmax = which(result$pi_gamma == max(result$pi_gamma))[1]
    result$effect = Paras0$gamma
    Paras0$gamma[argmax] = 0
    Paras0$paraUpdate$gamma[argmax] = FALSE

    # Fitting in the null hypothesis
    results0 = pleiClassify(Paras0)

    LR = result$LH-results0$LH
    LR[is.na(LR)] = 1
    chi2 = 2 * sum(LR)

    result$Pvalue = pchisq(chi2,1,lower.tail = F)
    result$effect = result$gamma[argmax]
    return(result)
}


# testing random effect of correlated pleiotropy
PCMR_testCorPlei = function(result){
    # corrected Pvalue range
    result$CHVP_test = pchisq(result$D_HVP,result$c,lower.tail = F)
    result$CHVP_test_Range = pchisq(result$D_HVP,
                                    result$c + c(-1.96,1.96) * result$c_sd,
                                    lower.tail = F)
    return(result)
}




