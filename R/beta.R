#' @rdname beta
#' @export
beta <- function(x,...){ UseMethod("beta",x) }
#' @rdname beta
#' @export
beta.default <- function(x,...){ stop('No default method.') }
#' @rdname beta
#' @export
beta.simplex <- function(x,num,den,a,...){
    snames <- names(x)
    out <- list()
    for (sname in snames){
        out[[sname]] <- beta(x=x[[sname]],num=num,den=den,a=a[[sname]],...)
    }
    out
}
#' @rdname alpha
#' @export
beta.standards <- function(x,num=num,den=den,a=a,...){
    beta(x=x$x,num=num,den=den,a=a,...)
}
#' @rdname beta
#' @export
beta.spot <- function(x,num,den,a,plot=FALSE,...){
    B <- common_denominator(c(num,den))
    groups <- groupbypairs(B)
    init <- init_beta(spot=x,groups=groups,a=a)
    fit <- optim(par=init,f=SS_b0g,method='L-BFGS-B',
                 lower=init-1,upper=init+1,spot=x,
                 groups=groups,a=a,hessian=TRUE)
    # TODO: covariance matrix
    if (plot){
        plot_beta(spot=x,groups=groups,b0g=out,a=a,...)
    }
    out
}

common_denominator <- function(ions){
    count <- table(ions)
    i <- which.max(count)
    out <- list()
    out$den <- names(count)[i]
    out$num <- names(count)[-i]
    out
}

init_beta <- function(spot,groups,a){
    b0 <- NULL
    g <- NULL
    b0names <- NULL
    gnames <- NULL
    for (gr in groups){
        b0 <- c(b0,log(a['exp_a0',gr$num]/a['exp_a0',gr$den]))
        b0names <- c(b0names,paste0('b0[',gr$num,'/',gr$den,']'))
        if (gr$elements['num']!=gr$elements['den']){
            g <- c(g,0)
            gnames <- c(gnames,paste0('g[',gr$elements['num'],
                                      '/',gr$elements['den'],']'))
        }
    }
    out <- c(b0,g)
    names(out) <- c(b0names,gnames)
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
        X <- Dp$t
        Y <- (Np$sig-Np$bkg)/(Dp$sig-Dp$bkg)
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

SS_b0g <- function(b0g,spot,groups,a){
    out <- 0
    for (gr in groups){
        bnames <- paste0('b0[',gr$num,'/',gr$den,']')
        b0 <- b0g[bnames]
        if (gr$elements['num']==gr$elements['den']){
            g <- 0
        } else {
            gname <- paste0('g[',gr$elements['num'],'/',gr$elements['den'],']')
            g <- b0g[gname]
        }
        ni <- length(gr$num)
        for (i in 1:ni){
            Np <- betapars(spot=spot,ion=gr$num[i],a=a)
            Dp <- betapars(spot=spot,ion=gr$den[i],a=a)
            pND <- predict_ND(b0=b0[i],g=g,Np=Np,Dp=Dp)
            SS <- sum((Np$sig - pND$N - Np$bkg)^2 +
                      (Dp$sig - pND$D - Dp$bkg)^2)
            out <- out + SS
        }
    }
    out
}

predict_ND <- function(b0,g,Np,Dp){
    out <- list()
    b0i <- b0 + g*Dp$t + Np$g*(Np$t-Dp$t)
    out$D <- get_exp_a0D(b0=b0,g=g,Np=Np,Dp=Dp)
    out$N <- out$D*exp(b0i)
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

# B = output of common_denominator
groupbypairs <- function(B){
    out <- list()
    out$den <- B$den
    out$num <- list()
    for (ion in B$num){
        nele <- element(ion)
        if (nele %in% names(out$num)){
            out$num[[nele]] <- c(out$num[[nele]],ion)
        } else {
            out$num[[nele]] <- ion
        }
    }
    out
}
