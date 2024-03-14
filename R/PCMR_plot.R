
PCMR_plot = function(result,Color = NA){

    # set color
    if(is.na(Color)[1]){
        Color = c("#4C5FFF","#949494") # [color_cause,color_correlated]
    }
    ColorInterval = adjustcolor(Color,alpha.f=0.35)

    # identify classify of IVs
    PleiClass = rep(0,length(result$LH))
    for(type in seq(result$Paras$num_gamma)){
        PleiClass = PleiClass + type * rowSums(result$W_ind[,type,])
    }

    # identify which is causal effect
    arg = which(result$pi_gamma == max(result$pi_gamma))[1]

    # the color of causal effect always be first index of Color.
    if(arg != 1){
        t = Color[arg]
        Color[arg] = Color[1]
        Color[1] = t

        t = ColorInterval[arg]
        ColorInterval[arg] = ColorInterval[1]
        ColorInterval[1] = t
    }

    # plot IVs by beta_exposure and beta_outcome
    xmin = min(result$Paras$beta_ex - 2 * result$Paras$beta_ex_se)
    xmax = max(result$Paras$beta_ex + 2 * result$Paras$beta_ex_se)
    ymin = min(result$Paras$beta_out - 2 * result$Paras$beta_out_se)
    ymax = max(result$Paras$beta_out + 2 * result$Paras$beta_out_se)

    plot( result$Paras$beta_ex,result$Paras$beta_out,
          xlim=c(xmin,xmax),
          ylim=c(ymin,ymax),
          xlab="Genetic association with exposure X",ylab="Genetic association with outcome Y",
          col=Color[PleiClass],font.lab = 3,pch=16)

    # plot conditional intervals of IVs
    for(i in seq(length(result$Paras$beta_ex))){
        lines(c(result$Paras$beta_ex[i] - 1.96 * result$Paras$beta_ex_se[i],result$Paras$beta_ex[i] + 1.96 * result$Paras$beta_ex_se[i]),c(result$Paras$beta_out[i],result$Paras$beta_out[i]),col=Color[PleiClass[i]])
        lines(c(result$Paras$beta_ex[i],result$Paras$beta_ex[i]),c(result$Paras$beta_out[i] - 1.96 * result$Paras$beta_out_se[i],result$Paras$beta_out[i] + 1.96 * result$Paras$beta_out_se[i]),col=Color[PleiClass[i]])
    }

    # x-axis and y-axis
    abline(h=0,v=0)

    # plot the fit line of causal effect and the sum of causal and correlated pleiotropy effect


    for(i in seq(length(result$gamma))){
        abline(a=result$b,b=result$gamma[i],lwd=2,col=Color[i])

        x = c(seq(xmin* 1.1,xmax* 1.1,0.0001),
              seq(xmax* 1.1,xmin* 1.1,-0.0001))

        y1 = x[1:(length(x)%/%2)] * (result$gamma[i] - 1.96*result$gSigma2[i] ** 0.5) + result$b
        y2 = x[(length(x)%/%2+1):length(x)]* (result$gamma[i] + 1.96*result$gSigma2[i] ** 0.5) + result$b
        y = c(y1,y2)

        x[x > xmax * 1.1] = xmax* 1.1
        x[x < xmin* 1.1] = xmin* 1.1
        y[y > ymax* 1.1] = ymax* 1.1
        y[y < ymin* 1.1] = ymin* 1.1
        polygon(x,y, col=ColorInterval[i] ,border=FALSE)
    }

}
