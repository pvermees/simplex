#' @rdname beta
#' @export
beta <- function(x,...){ UseMethod("beta",x) }
#' @rdname beta
#' @export
beta.default <- function(x,...){ stop('No default method.') }
#' @rdname beta
#' @export
beta.spot <- function(x,num,den,a,plot=FALSE,...){
    groups <- groupbypairs(num,den)
    out <- NULL
    if (plot){
        plot_beta(spot=x,num=num,den=den,b0g=NULL,a=a,...)
    }
    for (gr in groups){
        ni <- length(gr$num)
        init <- a['a0',gr$num]-a['a0',gr$den]
        if (gr$elements['num']==gr$elements['den']){
            fit <- optim(par=init,f=SS_b0,method='L-BFGS-B',
                         lower=init-1,upper=init+1,spot=x,
                         num=gr$num,den=gr$den,a=a)
        } else {
            init <- c(init,0)
            fit <- optim(par=init,f=SS_b0g,method='L-BFGS-B',
                         lower=init-1,upper=init+1,spot=x,
                         num=gr$num,den=gr$den,a=a)
        }
        out <- cbind(out)
    }
    rownames(out) <- c('b0','g')
    colnames(out) <- paste(num,den,sep='/')
    out
}

plot_beta <- function(spot,num,den,b0g,a,...){
    np <- length(num) # number of plot panels
    nr <- ceiling(sqrt(np)) # number of rows
    nc <- ceiling(np/nr)    # number of columns
    oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    for (i in 1:np){
        Np <- betapars(spot=spot,ion=num[i],a=a)
        Dp <- betapars(spot=spot,ion=den[i],a=a)
        driftcor <- exp(Np$g*(Np$t-Dp$t))
        X <- Dp$t
        Y <- driftcor*(Np$sig-Np$bkg)/(Dp$sig-Dp$bkg)
        ylab <- paste0('(',num[i],'-b)/(',den[i],'-b)')
        graphics::plot(X,Y,type='p',xlab='',ylab='',...)
        graphics::mtext(side=1,text='t',line=2)
        graphics::mtext(side=2,text=ylab,line=2)
    }
    graphics::par(oldpar)
}

SS_b0g <- function(b0g,spot,num,den,a){
    print(b0g)
    ni <- length(num)
    nt <- nrow(spot$time)
    out <- 0
    for (i in 1:ni){
        a0 <- rep(a['a0',num],nt) # TODO
        fit <- optim(a0,SS_a0,method='BFGS',b0g=b0g,spot=spot,
                     num=num[i],den=den[i],a=a)
        out <- out + fit$value
    }
    out
}
SS_b0 <- function(b0,spot,num,den,a){
    SS_b0g(c(b0,0),spot=spot,num=num,den=den,a=a)
}

SS_a0 <- function(a0,b0g,spot,num,den,a){
    nb0 <- length(b0g)-1
    b0 <- b0g[1:nb0]
    g <- b0g[nb0+1]
    Np <- betapars(spot=spot,ion=num,a=a)
    Dp <- betapars(spot=spot,ion=den,a=a)
    a0D <- a0
    a0N <- a0D + b0 + g*Dp$t + Np$sig*(Np$t-Dp$t)
    SS <- (Np$sig-Np$bkg-exp(a0N))^2 + (Dp$sig-Dp$bkg-exp(a0D))^2
    sum(SS)
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
            alreadyin <- (nele[i] %in% out[1,] || dele[i] %in% out[2,])
            if (!alreadyin){
                out <- cbind(out,c(nele[i],dele[i]))
            }
        }
    }
    out
}
