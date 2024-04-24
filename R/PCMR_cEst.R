PCMR_cEst = function(result,ref_beta_outcome,ref_se_outcome,samples=100,sample_boot = 30,cores=15, err_boot=100){


    result$est1 = PCMR(result$Paras$beta_ex, result$Paras$beta_ex_se,
                       result$Paras$beta_out,result$Paras$beta_out_se,num_gamma = 1,model="1")$gamma

    # Prepare for multithread processing
    Seqs = list()
    for(i in seq(samples*sample_boot)){
        if(i %% sample_boot == 1){
            Seqs[[i]] = sample(seq(length(ref_beta_outcome)),(length(result$Paras$beta_ex)))
        }else{
            Seqs[[i]] = Seqs[[i-1]]
        }
    }

    # bootstrapping
    # multithread
    sub_cl <- makeCluster(cores)
    clusterExport(sub_cl,"Seqs",envir = environment())
    clusterExport(sub_cl,".boot",envir = environment())
    Cup = clusterApplyLB(sub_cl,seq(samples*sample_boot),.boot_sample,result,ref_beta_outcome,ref_se_outcome)
    Cup_gamma = matrix(unlist(Cup),ncol=2,byrow=T)
    stopCluster(sub_cl) # close

    # Data decomposition
    Cup_chi2 = c()
    for(i in seq(samples)){
        Cup_gamma_sub = Cup_gamma[seq((i-1)*sample_boot+1,i*sample_boot),]
        Cup_chi2 = c(Cup_chi2,(mean(Cup_gamma_sub[,1]) - mean(Cup_gamma_sub[,2]))^2/( var(Cup_gamma_sub[,1]) +  var(Cup_gamma_sub[,2])))
    }

    # fitting to estimate parameter c, and the range error
    temp = function(c,Chi2s){
        return(sum((sort(-log(pchisq(Cup_chi2,1*c,lower.tail = F))) - sort(-log(seq(Cup_chi2)/length(Cup_chi2))))^2))
    }

    result$c = optimize(temp,lower=0,upper=10,Chi2s=Cup_chi2)$minimum

    # estimate standard error of estimate correct_factor
    C_tmp = c()
    for(i in seq(err_boot)){
        C_tmp = c(C_tmp,optimize(temp,lower=0,upper=10,Chi2s=sample(Cup_chi2,replace = T))$minimum)
    }
    result$c_sd = sd(C_tmp)

    return(result)
}
