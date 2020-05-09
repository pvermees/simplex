#' @rdname beta
#' @export
beta <- function(x,...){ UseMethod("beta",x) }
#' @rdname beta
#' @export
beta.default <- function(x,...){ stop('No default method.') }
#' @rdname beta
#' @export
beta.spot <- function(x,num,den,a,plot=FALSE,...){
    out <- rbind(log(a['exp_a0',num]/a['exp_a0',den]),0) # initialise
    rownames(out) <- c('b0','g')
    colnames(out) <- paste0(num,'/',den)
    groups <- groupbypairs(num,den)
    for (gr in groups){
        print(gr$num)
        nb <- length(gr$num)
        ratios <- paste0(gr$num,'/',gr$den)
        init <- out['b0',ratios]
        if (gr$elements['num']==gr$elements['den']){
            fit <- optim(par=init,f=SS_b0,method='L-BFGS-B',
                         lower=init-1,upper=init+1,spot=x,
                         num=gr$num,den=gr$den,a=a)
        } else {
            init <- c(init,0)
            fit <- optim(par=init,f=SS_b0g,method='L-BFGS-B',
                         lower=init-1,upper=init+1,spot=x,
                         num=gr$num,den=gr$den,a=a)
            out['g',ratios] <- fit$par[nb+1]
        }
        out['b0',ratios] <- fit$par[1:nb]
    }
    if (plot){
        plot_beta(spot=x,num=num,den=den,b0g=out,a=a,...)
    }
    out
}

plot_beta <- function(spot,num,den,b0g,a,...){
    np <- length(num) # number of plot panels
    nr <- ceiling(sqrt(np)) # number of rows
    nc <- ceiling(np/nr)    # number of columns
    nt <- nrow(spot$time)
    oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    for (i in 1:np){
        ratio <- paste0(num[i],'/',den[i])
        Np <- betapars(spot=spot,ion=num[i],a=a)
        Dp <- betapars(spot=spot,ion=den[i],a=a)
        driftcor <- exp(Np$g*(Np$t-Dp$t))
        X <- Dp$t
        Y <- driftcor*(Np$sig-Np$bkg)/(Dp$sig-Dp$bkg)
        b0 <- b0g['b0',ratio]
        g <- b0g['g',ratio]
        pND <- predict_ND(b0=b0,g=g,Np=Np,Dp=Dp)
        Ypred <- pND$N/pND$D
        ylab <- paste0('(',num[i],'-b)/(',den[i],'-b)')
        graphics::plot(c(X,X),c(Y,Ypred),type='n',xlab='',ylab='',...)
        graphics::points(X,Y)
        graphics::lines(X,Ypred)
        graphics::mtext(side=1,text='t',line=2)
        graphics::mtext(side=2,text=ylab,line=2)
    }
    graphics::par(oldpar)
}

SS_b0g <- function(b0g,spot,num,den,a){
    ni <- length(num)
    b0 <- b0g[1:ni]
    g <- b0g[ni+1]
    out <- 0
    for (i in 1:ni){
        Np <- betapars(spot=spot,ion=num[i],a=a)
        Dp <- betapars(spot=spot,ion=den[i],a=a)
        pND <- predict_ND(b0=b0[i],g=g,Np=Np,Dp=Dp)
        SS <- sum((Np$sig - pND$N)^2 + (Dp$sig - pND$D)^2)
        out <- out + SS
    }
    out
}

predict_ND <- function(b0,g,Np,Dp){
    exp_a0D <- get_exp_a0D(b0=b0,g=g,Np=Np,Dp=Dp)
    b0i <- b0 + g*Dp$t + Np$g*(Np$t-Dp$t)
    exp_a0N <- exp_a0D*exp(b0i)
    out <- list()
    out$N <- Np$bkg + exp_a0N
    out$D <- Dp$bkg + exp_a0D
    out
}

SS_b0 <- function(b0,spot,num,den,a){
    SS_b0g(c(b0,0),spot=spot,num=num,den=den,a=a)
}

get_exp_a0D <- function(b0,g,Np,Dp){
    b0i <- b0 + g*Dp$t + Np$g*(Np$t-Dp$t)
    numerator <- (Dp$sig-Dp$bkg) + (Np$sig-Np$bkg)*exp(b0i)
    denominator <- exp(b0i)^2 + 1
    numerator/denominator
}

betapars <- function(spot,ion,a){
    out <- list()
    out$sig <- spot$signal[,ion]
    out$bkg <- background(spot,ion)
    out$g <- a['g',ion]
    out$t <- hours(spot$time[,ion])
    out
}

groupbypairs <- function(num,den){
    out <- list()
    pairs <- elementpairs(num,den)
    nele <- element(num)
    dele <- element(den)
    for (i in 1:ncol(pairs)){
        matches <- which((nele %in% pairs[1,i]) & (dele %in% pairs[2,i]))
        out[[i]] <- list()
        out[[i]]$elements <- pairs[,i]
        names(out[[i]]$elements) <- c('num','den')
        out[[i]]$num <- num[matches]
        out[[i]]$den <- den[matches]
    }
    out
}

elementpairs <- function(num,den){
    nele <- element(num)
    dele <- element(den)
    out <- rbind(nele[1],dele[1])
    if (length(num)>1){
        for (i in 2:length(num)){
            ni <- which(out[1,] %in% nele[i])
            di <- which(out[2,] %in% dele[i])
            alreadyin <- any(ni%in%di)
            if (!alreadyin){
                out <- cbind(out,c(nele[i],dele[i]))
            }
        }
    }
    out
}
