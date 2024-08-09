PCMR_testCausal_bootstrap = function(result,samples=1000,cores=15){

    # bootstrapping for observed IVs
    sub_cl <- makeCluster(cores)
    clusterExport(sub_cl,"result",envir = environment())

    Cup = parSapply(sub_cl,seq(samples),.boot,result,simplify = TRUE)
    stopCluster(sub_cl) # close

    Ncol = result$Paras$num_gamma + 1
    Cup_gamma = matrix(unlist(Cup),ncol=Ncol,byrow=T)

    gamma_min = Cup_gamma[,1]
    gamma_max = Cup_gamma[,Ncol - 1]

    result$bootstrap = list()

    result$bootstrap$mean_minClass = mean(gamma_min)
    result$bootstrap$mean_maxClass = mean(gamma_max)

    result$bootstrap$sd_minClass = sd(gamma_min)
    result$bootstrap$sd_maxClass = sd(gamma_max)

    result$bootstrap$min_interval = stringr::str_c(quantile(gamma_min,c(0.025,0.975)),collapse = "_")
    result$bootstrap$max_interval = stringr::str_c(quantile(gamma_max,c(0.025,0.975)),collapse = "_")

    result$bootstrap$P_min = pt(abs(mean(gamma_min)/result$bootstrap$sd_minClass), df=length(result$Paras$beta_ex)-1, lower.tail=F)*2
    result$bootstrap$P_max = pt(abs(mean(gamma_max)/result$bootstrap$sd_maxClass), df=length(result$Paras$beta_ex)-1, lower.tail=F)*2

    result$bootstrap$P_min_permut = 1 - abs(sum(gamma_min > 0) - sum(gamma_min < 0))/samples
    result$bootstrap$P_max_permut = 1 - abs(sum(gamma_max > 0) - sum(gamma_max < 0))/samples

    result$D_HVP = (result$bootstrap$mean_minClass - result$bootstrap$mean_maxClass)^2/
                        ( result$bootstrap$sd_minClass^2 +  result$bootstrap$sd_maxClass^2)

    # The probability of dominant IV group
    result$dom_prob = c(table(Cup_gamma[,Ncol]))/sum(Cup_gamma[,Ncol]>0)
    result$dom_prob = result$dom_prob[order(result$gamma)]

    return(result)
}
