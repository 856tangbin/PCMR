PCMR_testCausal_bootstrap = function(result,samples=1000,cores=15){

    # bootstrapping for observed IVs
    sub_cl <- makeCluster(cores)
    clusterExport(sub_cl,"result",envir = environment())

    Cup = parSapply(sub_cl,seq(samples),.boot,result,simplify = TRUE)
    stopCluster(sub_cl) # close

    Cup_gamma = matrix(unlist(Cup),ncol=2,byrow=T)

    result$bootstrap = list()

    result$bootstrap$mean_minClass = mean(Cup_gamma[,1])
    result$bootstrap$mean_maxClass = mean(Cup_gamma[,2])

    result$bootstrap$sd_minClass = sd(Cup_gamma[,1])
    result$bootstrap$sd_maxClass = sd(Cup_gamma[,2])

    result$bootstrap$min_interval = stringr::str_c(quantile(Cup_gamma[,1],c(0.025,0.975)),collapse = "_")
    result$bootstrap$max_interval = stringr::str_c(quantile(Cup_gamma[,2],c(0.025,0.975)),collapse = "_")

    result$bootstrap$P_min = pt(abs(mean(Cup_gamma[,1])/result$bootstrap$sd_minClass), df=length(result$Paras$beta_ex)-1, lower.tail=F)*2
    result$bootstrap$P_max = pt(abs(mean(Cup_gamma[,2])/result$bootstrap$sd_maxClass), df=length(result$Paras$beta_ex)-1, lower.tail=F)*2

    result$bootstrap$P_min_permut = 1 - abs(sum(Cup_gamma[,1] > 0) - sum(Cup_gamma[,1] < 0))/samples
    result$bootstrap$P_max_permut = 1 - abs(sum(Cup_gamma[,2] > 0) - sum(Cup_gamma[,2] < 0))/samples

    result$D_HVP = (result$bootstrap$mean_minClass - result$bootstrap$mean_maxClass)^2/
                        ( result$bootstrap$sd_minClass^2 +  result$bootstrap$sd_maxClass^2)
    return(result)
}
