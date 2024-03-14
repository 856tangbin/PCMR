pleiClassify = function(Paras){

    # whether using the intact model?
    if(Paras$isIntact){
        return(pleiClassifyIntact(Paras))
    }

    size = length(Paras$beta_ex)

    # Initialization
    # gamma
    if(is.na(Paras$gamma[1])){
        gamma = quantile(Paras$beta_out/Paras$beta_ex,seq(Paras$num_gamma)/(Paras$num_gamma+1))
    }else{
        gamma = Paras$gamma
    }

    # b
    if(is.na(Paras$b[1])){
        b = 0
    }else{
        b = Paras$b
    }

    # pi_gamma
    if(is.na(Paras$pi_gamma[1])){
        pi_gamma = seq(Paras$num_gamma)
        pi_gamma = pi_gamma / sum(pi_gamma)
    }else{
        pi_gamma = Paras$pi_gamma
        pi_gamma = pi_gamma / sum(pi_gamma)
    }

    # pi_b
    if(is.na(Paras$pi_b[1])){
        pi_b = seq(Paras$num_b)
        pi_b = pi_b / sum(pi_b)
    }else{
        pi_b = Paras$pi_b
        pi_b = pi_b / sum(pi_b)
    }

    # gSigma2
    if(is.na(Paras$gSigma2[1])){
        gSigma2 = rep(0,Paras$num_gamma)
    }else{
        gSigma2 = Paras$gSigma2
    }

    # max_gSigma2
    if(is.na(Paras$max_gSigma2[1])){
        max_gSigma2 = 1e8
    }else{
        max_gSigma2 = Paras$max_gSigma2
    }

    # sigma2
    if(is.na(Paras$sigma2[1])){
        sigma2 = rep(0,Paras$num_b)
    }else{
        sigma2 = Paras$sigma2
    }

    # max sigma2
    if(is.na(Paras$max_sigma2[1])){
        max_sigma2 = 1e8
    }else{
        max_sigma2 = Paras$max_sigma2
    }


    Cup_gamma = rep(-1000,5)

    # EM Algorithm
    for(. in seq(Paras$maxIteration)){
        gamma_b = gamma[order(pi_gamma)]
        b_b = b
        sigma2_b = sigma2[order(pi_b)]
        gSigma2_b = gSigma2[order(pi_gamma)]

        # E-step
        W_ = array(0,dim=c(size,Paras$num_gamma,Paras$num_b))
        W = W_

        # Q function
        Qfun = function(ex,exSE,out,outSE){
            S = array(0,dim=c(Paras$num_gamma,Paras$num_b))
            for(j in seq(Paras$num_gamma)){
                for(k in seq(Paras$num_b)){
                    S[j,k] = dnorm(out - gamma[j] * ex,mean=b,sd=sqrt(outSE**2 + ex**2 * gSigma2[j] +sigma2[k]))
                }
            }
            return(S)
        }

        for(i in seq(size)){
            Q = Qfun(Paras$beta_ex[i], Paras$beta_ex_se[i], Paras$beta_out[i], Paras$beta_out_se[i])
            W_[i,,] = Q
            if(sum(Q) != 0){
                W[i,,] = Q/sum(Q)
            }else{
                W[i,,] = array(0,dim=dim(Q))
            }

        }


        # M-step of gamma
        for(j in seq(Paras$num_gamma)){
            num = 0
            den = 0
            for(k in seq(Paras$num_b)){
                num = num + sum(W[,j,k] * Paras$beta_ex * (Paras$beta_out - b)/(Paras$beta_out_se^2 + Paras$beta_ex**2 * gSigma2[j] + sigma2[k]))
                den = den + sum(W[,j,k] * Paras$beta_ex ^ 2 / (Paras$beta_out_se ^ 2 + Paras$beta_ex**2 * gSigma2[j] + sigma2[k]))
            }
            if (Paras$paraUpdate$gamma[j]){
                if(den != 0){
                    gamma[j] = num/den
                }else{
                    gamma[j] = num/Paras$thred
                }
            }
        }

        # M-step of gSigma2
        for(j in seq(Paras$num_gamma)){
            f = function(gSigma2){
                t = 0
                for(k in seq(Paras$num_b)){
                    t = t + sum(W[,j,k]*log(dnorm(Paras$beta_out - gamma[j] * Paras$beta_ex,mean=b,sd=sqrt(Paras$beta_out_se**2 + Paras$beta_ex**2 * gSigma2 + sigma2[k])) * pi_gamma[j] * pi_b[k]/W[,j,k]))
                }

                return(-t)
            }
            if(Paras$paraUpdate$gSigma2[j]){
                gSigma2[j] = optimize(f,lower=-Paras$thred**2,upper=max_gSigma2,tol=Paras$thred**2)$minimum

                if(gSigma2[j] < 0){
                    gSigma2[j] = 0
                }
            }
        }

        # M-step of b
        num = 0
        den = 0
        for(k in seq(Paras$num_b)){
            for(j in seq(Paras$num_gamma)){
                num = num + sum(W[,j,k] * (Paras$beta_out - gamma[j] * Paras$beta_ex) / (Paras$beta_out_se ^ 2 + Paras$beta_ex**2 * gSigma2[j] + sigma2[k]))
                den = den + sum(W[, j, k] / (Paras$beta_out_se ^ 2 + Paras$beta_ex**2 * gSigma2[j] + sigma2[k]))
            }
        }

        if(Paras$paraUpdate$b[1]){
            if(den != 0){
                b = num/den
            }else{
                b = num/Paras$thred
            }
        }


        # M-step of pi_gamma
        for(j in seq(Paras$num_gamma)){
            if(Paras$paraUpdate$pi_gamma[j]){
                pi_gamma[j] = sum(W[,j,])/size
            }
        }

        # M-step of pi_b
        for(k in seq(Paras$num_b)){
            if(Paras$paraUpdate$pi_b[k]){
                pi_b[k] = sum(W[,,k])/size
            }
        }

        for(k in seq(Paras$num_b)){
            f = function(sigma2){
                t = 0
                for(j in seq(Paras$num_gamma)){
                    t = t + sum(W[,j,k]*log(dnorm(Paras$beta_out - gamma[j] * Paras$beta_ex,mean=b,sd=sqrt(Paras$beta_out_se**2 + Paras$beta_ex**2 * gSigma2[j] + sigma2)) * pi_gamma[j] * pi_b[k]/W[,j,k]))
                }
                return(-t)
            }
            if(Paras$paraUpdate$sigma2[k]){
                sigma2[k] = optimize(f,lower=-Paras$thred**2,upper=max_sigma2,tol=Paras$thred**2)$minimum

                if(sigma2[k] < 0){
                    sigma2[k] = 0
                }
            }
        }

        # Convergence judgment
        if(sum(abs(gamma_b - gamma[order(pi_gamma)])) < Paras$thred &
           sum(abs(b_b - b)) < Paras$thred &
           sum(abs(sigma2_b - sigma2[order(pi_b)])) < Paras$thred**2 &
           sum(abs(gSigma2_b - gSigma2[order(pi_gamma)])) < Paras$thred**2){
            break
        }

        # extreme convergence case
        if(gamma[1] %in% Cup_gamma){
            break
        }
        Cup_i = . %% length(Cup_gamma)
        Cup_gamma[Cup_i] = gamma[1]

    }

    W_ = array(0,dim=c(size,Paras$num_gamma,Paras$num_b))
    W = W_

    for(i in seq(size)){
        Q = Qfun(Paras$beta_ex[i], Paras$beta_ex_se[i], Paras$beta_out[i], Paras$beta_out_se[i])
        W_[i,,] = Q
        if(sum(Q) != 0){
            W[i,,] = Q/sum(Q)
        }else{
            W[i,,] = array(0,dim=dim(Q))
        }
    }

    # log-liklihood
    LH = rep(0,size)
    for(i in seq(size)){
        LH[i] = sum(log(t(t(matrix(W_[i,,],nrow = Paras$num_gamma,ncol = Paras$num_b)) * pi_b) * pi_gamma/W[i,,]) * W[i,,])
    }

    # AIC
    AIC = 2 * (numOfPara(Paras$paraUpdate) - 2) - 2 * sum(LH) # -2 represents the parameters of pi_gamma and pi_b

    # Classify
    W_ind = W_ * 0
    for(i in seq(size)){
        arg = which(max(W_[i,,]) == W_[i,,])[1]
        W_ind[i,(arg - 1)%%Paras$num_gamma + 1,(arg-1)%/%(Paras$num_gamma)+1] = 1
    }

    return(list(Paras=Paras,
                gamma=gamma,pi_gamma=pi_gamma,b=b,pi_b = pi_b,sigma2 = sigma2,gSigma2=gSigma2,
                LH=LH,AIC=AIC,W_ind=W_ind,W=W))
}


pleiClassifyIntact = function(Paras){


    size = length(Paras$beta_ex)
    min_gamma = min(Paras$beta_out/Paras$beta_ex)
    max_gamma = max(Paras$beta_out/Paras$beta_ex)


    # Initialization
    # gamma
    if(is.na(Paras$gamma[1])){
        gamma = quantile(Paras$beta_out/Paras$beta_ex,seq(Paras$num_gamma)/(Paras$num_gamma+1))
    }else{
        gamma = Paras$gamma
    }

    # b
    if(is.na(Paras$b[1])){
        b = 0
    }else{
        b = Paras$b
    }

    # pi_gamma
    if(is.na(Paras$pi_gamma[1])){
        pi_gamma = seq(Paras$num_gamma)
        pi_gamma = pi_gamma / sum(pi_gamma)
    }else{
        pi_gamma = Paras$pi_gamma
        pi_gamma = pi_gamma / sum(pi_gamma)
    }

    # pi_b
    if(is.na(Paras$pi_b[1])){
        pi_b = seq(Paras$num_b)
        pi_b = pi_b / sum(pi_b)
    }else{
        pi_b = Paras$pi_b
        pi_b = pi_b / sum(pi_b)
    }

    # gSigma2
    if(is.na(Paras$gSigma2[1])){
        gSigma2 = rep(0,Paras$num_gamma)
    }else{
        gSigma2 = Paras$gSigma2
    }

    # max_gSigma2
    if(is.na(Paras$max_gSigma2[1])){
        max_gSigma2 = 1e8
    }else{
        max_gSigma2 = Paras$max_gSigma2
    }

    # sigma2
    if(is.na(Paras$sigma2[1])){
        sigma2 = rep(0,Paras$num_b)
    }else{
        sigma2 = Paras$sigma2
    }

    # max sigma2
    if(is.na(Paras$max_sigma2[1])){
        max_sigma2 = 1e8
    }else{
        max_sigma2 = Paras$max_sigma2
    }




    # EM Algorithm
    for(. in seq(Paras$maxIteration)){
        gamma_b = gamma[order(pi_gamma)]
        b_b = b
        sigma2_b = sigma2[order(pi_b)]
        gSigma2_b = gSigma2[order(pi_gamma)]

        # E-step
        W_ = array(0,dim=c(size,Paras$num_gamma,Paras$num_b))
        W = W_

        # Q function
        Qfun = function(ex,exSE,out,outSE){
            S = array(0,dim=c(Paras$num_gamma,Paras$num_b))
            for(j in seq(Paras$num_gamma)){
                for(k in seq(Paras$num_b)){
                    S[j,k] = dnorm(out - gamma[j] * ex,mean=b,sd=sqrt(outSE**2 +
                                                                          exSE**2 * gamma[j]**2 - 2 * gamma[j] * Paras$rho * exSE * outSE
                                                                      + ex**2 * gSigma2[j] + sigma2[k]))
                }
            }
            return(S)
        }

        for(i in seq(size)){
            Q = Qfun(Paras$beta_ex[i], Paras$beta_ex_se[i], Paras$beta_out[i], Paras$beta_out_se[i])
            W_[i,,] = Q
            if(sum(Q) != 0){
                W[i,,] = Q/sum(Q)
            }else{
                W[i,,] = array(0,dim=dim(Q))
            }

        }


        # M-step of gamma
        for(j in seq(Paras$num_gamma)){
            f = function(gamma){
                t = 0
                for(k in seq(Paras$num_b)){
                    t = t + sum(W[,j,k]*log((dnorm(Paras$beta_out - gamma * Paras$beta_ex,mean=b,sd=sqrt(Paras$beta_out_se**2 +
                                                                                                 Paras$beta_ex_se**2 * gamma**2 - 2 * gamma * Paras$rho * Paras$beta_ex_se * Paras$beta_out_se
                                                                                             + Paras$beta_ex**2 * gSigma2[j] + sigma2[k]))) * pi_gamma[j] * pi_b[k]/W[,j,k]))
                }

                return(-t)
            }
            if(Paras$paraUpdate$gamma[j]){
                gamma[j] = optimize(f,lower=min_gamma,upper=max_gamma,tol=Paras$thred**2)$minimum
            }
        }

        # M-step of gSigma2
        for(j in seq(Paras$num_gamma)){
            f = function(gSigma2){
                t = 0
                for(k in seq(Paras$num_b)){
                    t = t + sum(W[,j,k]*log(dnorm(Paras$beta_out - gamma[j] * Paras$beta_ex,mean=b,sd=sqrt(Paras$beta_out_se**2 +
                                                                                                   Paras$beta_ex_se**2 * gamma[j]**2 - 2 * gamma[j] * Paras$rho * Paras$beta_ex_se * Paras$beta_out_se
                                                                                               + Paras$beta_ex**2 * gSigma2 + sigma2[k])) * pi_gamma[j] * pi_b[k]/W[,j,k]))
                }

                return(-t)
            }
            if(Paras$paraUpdate$gSigma2[j]){
                gSigma2[j] = optimize(f,lower=-Paras$thred**2,upper=max_gSigma2,tol=Paras$thred**2)$minimum

                if(gSigma2[j] < 0){
                    gSigma2[j] = 0
                }
            }
        }

        # M-step of b
        num = 0
        den = 0
        for(k in seq(Paras$num_b)){
            for(j in seq(Paras$num_gamma)){
                num = num + sum(W[,j,k] * (Paras$beta_out - gamma[j] * Paras$beta_ex) / (Paras$beta_out_se ^ 2 +
                                                                                 Paras$beta_ex_se**2 * gamma[j]**2 - 2 * gamma[j] * Paras$rho * Paras$beta_ex_se * Paras$beta_out_se
                                                                             + Paras$beta_ex**2 * gSigma2[j] + sigma2[k]))
                den = den + sum(W[, j, k] / (Paras$beta_out_se ^ 2 +
                                                 Paras$beta_ex_se**2 * gamma[j]**2 - 2 * gamma[j] * Paras$rho * Paras$beta_ex_se * Paras$beta_out_se
                                             + Paras$beta_ex**2 * gSigma2[j] + sigma2[k]))
            }
        }

        if(Paras$paraUpdate$b[1]){
            if(den != 0){
                b = num/den
            }else{
                b = num/Paras$thred
            }
        }


        # M-step of pi_gamma
        for(j in seq(Paras$num_gamma)){
            if(Paras$paraUpdate$pi_gamma[j]){
                pi_gamma[j] = sum(W[,j,])/size
            }
        }

        # M-step of pi_b
        for(k in seq(Paras$num_b)){
            if(Paras$paraUpdate$pi_b[k]){
                pi_b[k] = sum(W[,,k])/size
            }
        }

        for(k in seq(Paras$num_b)){
            f = function(sigma2){
                t = 0
                for(j in seq(Paras$num_gamma)){
                    t = t + sum(W[,j,k]*log(dnorm(Paras$beta_out - gamma[j] * Paras$beta_ex,mean=b,sd=sqrt(Paras$beta_out_se**2 +
                                                                                                   Paras$beta_ex_se**2 * gamma[j]**2 - 2 * gamma[j] * Paras$rho * Paras$beta_ex_se * Paras$beta_out_se
                                                                                               + Paras$beta_ex**2 * gSigma2[j] + sigma2)) * pi_gamma[j] * pi_b[k]/W[,j,k]))
                }
                return(-t)
            }
            if(Paras$paraUpdate$sigma2[k]){
                sigma2[k] = optimize(f,lower=-Paras$thred**2,upper=max_sigma2,tol=Paras$thred**2)$minimum

                if(sigma2[k] < 0){
                    sigma2[k] = 0
                }
            }
        }

        # Convergence judgment
        if(sum(abs(gamma_b - gamma[order(pi_gamma)])) < Paras$thred &
           sum(abs(b_b - b)) < Paras$thred &
           sum(abs(sigma2_b - sigma2[order(pi_b)])) < Paras$thred**2 &
           sum(abs(gSigma2_b - gSigma2[order(pi_gamma)])) < Paras$thred**2){
            break
        }
    }

    W_ = array(0,dim=c(size,Paras$num_gamma,Paras$num_b))
    W = W_

    for(i in seq(size)){
        Q = Qfun(Paras$beta_ex[i], Paras$beta_ex_se[i], Paras$beta_out[i], Paras$beta_out_se[i])
        W_[i,,] = Q
        if(sum(Q) != 0){
            W[i,,] = Q/sum(Q)
        }else{
            W[i,,] = array(0,dim=dim(Q))
        }
    }

    # log-liklihood
    LH = rep(0,size)
    for(i in seq(size)){
        LH[i] = sum(log(t(t(matrix(W_[i,,],nrow = Paras$num_gamma,ncol = Paras$num_b)) * pi_b) * pi_gamma/W[i,,]) * W[i,,])
    }

    # AIC
    AIC = 2 * (numOfPara(Paras$paraUpdate) - 2) - 2 * sum(LH) # -2 represents the parameters of pi_gamma and pi_b

    # Classify
    W_ind = W_ * 0
    for(i in seq(size)){
        arg = which(max(W_[i,,]) == W_[i,,])[1]
        W_ind[i,(arg - 1)%%Paras$num_gamma + 1,(arg-1)%/%(Paras$num_gamma)+1] = 1
    }

    return(list(Paras=Paras,
                gamma=gamma,pi_gamma=pi_gamma,b=b,pi_b = pi_b,sigma2 = sigma2,gSigma2=gSigma2,
                LH=LH,AIC=AIC,W_ind=W_ind,W=W))
}

numOfPara = function(paraUpdate){
    S = 0
    for(name in names(paraUpdate)){
        S = S + sum(paraUpdate[[name]])
    }
    return(S)
}

