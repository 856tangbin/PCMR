.boot = function(.,result){
    library(PCMR)

    set.seed(randomSeeds[.])

    size = length(result$Paras$beta_ex)
    Ind = sample(seq(size),size = size,replace = T)

    Paras = result$Paras
    Paras$beta_ex = result$Paras$beta_ex[Ind]
    Paras$beta_ex_se = result$Paras$beta_ex_se[Ind]
    Paras$beta_out = result$Paras$beta_out[Ind]
    Paras$beta_out_se = result$Paras$beta_out_se[Ind]
    results0 = pleiClassify(Paras)


    if (sd(results0$pi_gamma) < results0$Paras$thred){
        dominant_group = -1
    }else{
        largest_group = order(results0$pi_gamma,decreasing = T)[1]
        for(ind in order(results0$gamma)){
            if(ind == largest_group){
                dominant_group = ind
            }
        }
    }

    return(c(sort(results0$gamma),dominant_group))
}

.boot_sample = function(.,result,beta_outcome,se_outcome){
    library(PCMR)

    set.seed(randomSeeds[.])

    result_sample = result
    result_sample$Paras$beta_out = beta_outcome[Seqs[[.]]] + result$est1 * rnorm(length(result$Paras$beta_ex),result$Paras$beta_ex,result$Paras$beta_ex_se)
    result_sample$Paras$beta_out_se = se_outcome[Seqs[[.]]]

    return(.boot(.,result_sample))
}
