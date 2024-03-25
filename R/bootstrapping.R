.boot = function(.,result){
    library(PCMR)

    size = length(result$Paras$beta_ex)
    Ind = sample(seq(size),size = size,replace = T)

    Paras = result$Paras
    Paras$beta_ex = result$Paras$beta_ex[Ind]
    Paras$beta_ex_se = result$Paras$beta_ex_se[Ind]
    Paras$beta_out = result$Paras$beta_out[Ind]
    Paras$beta_out_se = result$Paras$beta_out_se[Ind]
    results0 = pleiClassify(Paras)

    return(sort(results0$gamma))
}

.boot_sample = function(.,result,beta_outcome,se_outcome){
    library(PCMR)

    result_sample = result

    result_sample$Paras$beta_out = beta_outcome[Seqs[[.]]] + result$est1 * result$Paras$beta_ex
    result_sample$Paras$beta_out_se = sqrt(se_outcome[Seqs[[.]]]**2 + (result$est1 * result$Paras$beta_ex_se)**2)

    return(.boot(.,result_sample))
}
