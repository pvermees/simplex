#' @rdname logratios
#' @export
logratios <- function(x,...){ UseMethod("logratios",x) }
#' @rdname logratios
#' @export
logratios.default <- function(x,...){ stop('No default method.') }
#' @rdname logratios
#' @export
logratios.simplex <- function(x,num,den,dc=NULL,...){
    snames <- names(x)
    out <- list()
    for (sname in snames){
        print(sname)
        out[[sname]] <- logratios(x=x[[sname]],num=num,den=den,dc=dc[[sname]],...)
    }
    out
}
#' @rdname alpha
#' @export
logratios.standards <- function(x,num,den,dc=NULL,...){
    logratios(x=x$x,num=num,den=den,dc=dc,...)
}
#' @rdname logratios
#' @export
logratios.spot <- function(x,num,den,dc=NULL,plot=FALSE,...){
    B <- common_denominator(c(num,den))
    groups <- groupbypairs(B)
    init <- init_logratios(spot=x,groups=groups,dc=dc)
    fit <- optim(par=init,f=SS_b0g,method='L-BFGS-B',
                 lower=init-2,upper=init+2,spot=x,
                 groups=groups,dc=dc,hessian=TRUE)
    fit$mse <- SS_b0g(fit$par,spot=x,groups=groups,dc=dc,mse=TRUE)
    out <- common2original(fit=fit,num=num,den=den,groups=groups)
    if (plot){
        plot_logratios(spot=x,groups=groups,b0g=out$b0g,dc=dc,...)
    }
    out
}

# converts logratio intercepts and slopes from common
# denominator to scientifically useful logratios
common2original <- function(fit,num,den,groups){
    B0G <- b0g2list(b0g=fit$par,groups=groups)
    b0in <- B0G$b0
    gin <- B0G$g
    bnamesin <- names(b0in)
    gnamesin <- names(gin)
    rnames <- paste0(num,'/',den)
    outnames <- c(paste0('b0[',rnames,']'),
                  paste0('g[',rnames,']'))
    ni <- length(rnames)
    J <- matrix(0,nrow=2*ni,ncol=length(fit$par))
    colnames(J) <- names(fit$par)
    rownames(J) <- outnames
    b0gout <- rep(0,2*ni)
    names(b0gout) <- outnames
    for (i in 1:ni){
        nion <- num[i]
        dion <- den[i]
        if (nion == groups$den){
            iden <- which(bnamesin %in% dion)
            b0gout[i] <- b0gout[i] - b0in[iden]
            J[i,iden] <- -1
        } else {
            inum <- which(bnamesin %in% nion)
            b0gout[i] <- b0gout[i] + b0in[inum]
            J[i,inum] <- 1
            if (dion != groups$den){
                iden <- which(bnamesin %in% dion)
                b0gout[i] <- b0gout[i] - b0in[iden]
                J[i,iden] <- -1
            }
        }
        nele <- element(num[i])
        dele <- element(den[i])
        if (nele == dele){
            b0gout[ni+i] <- 0
        } else {
            inele <- which(gnamesin %in% nele)
            b0gout[ni+i] <- b0gout[ni+i] + gin[inele]
            J[ni+i,inele] <- 1
            if (dele != element(groups$den)){
                idele <- which(gnamesin %in% dele)
                b0gout[ni+i] <- b0gout[ni+i] - gin[idele]
                J[ni+i,idele] <- -1
            }
        }
    }
    E <- MASS::ginv(fit$hessian/2)*fit$mse
    out <- list()
    out$b0g <- b0gout
    out$cov <- J %*% E %*% t(J)
    out
}

# uses the most used ion as a common denominator
common_denominator <- function(ions){
    count <- table(ions)
    i <- which.max(count)
    out <- list()
    out$num <- names(count)[-i]
    out$den <- names(count)[i]
    out
}

init_logratios <- function(spot,groups,dc=NULL){
    b0 <- NULL
    g <- NULL
    b0names <- NULL
    gnames <- NULL
    den <- groups$den
    for (num in groups$num){
        b0 <- c(b0,log(dc['exp_a0',num]/dc['exp_a0',den]))
        b0names <- c(b0names,num)
        nele <- unique(element(num))
        dele <- element(den)
        if (nele!=dele){
            g <- c(g,0)
            gnames <- c(gnames,nele)
        }
    }
    out <- c(b0,g)
    names(out) <- c(b0names,gnames)
    out
}

plot_logratios <- function(spot,groups,b0g,dc=NULL,...){
    np <- length(num) # number of plot panels
    nr <- ceiling(sqrt(np)) # number of rows
    nc <- ceiling(np/nr)    # number of columns
    nt <- nrow(spot$time)
    nb <- length(b0g)/2
    oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    for (i in 1:np){
        ratio <- paste0(num[i],'/',den[i])
        Np <- betapars(spot=spot,ion=num[i],dc=dc)
        Dp <- betapars(spot=spot,ion=den[i],dc=dc)
        X <- Dp$t
        Y <- (Np$sig-Np$bkg)/(Dp$sig-Dp$bkg)
        b0 <- b0g[paste0('b0[',ratio,']')]
        g <- b0g[paste0('g[',ratio,']')]
        Ypred <- exp(b0 + g*Dp$t + dc['g',num[i]]*(Np$t-Dp$t))
        ylab <- paste0('(',num[i],'-b)/(',den[i],'-b)')
        graphics::plot(c(X,X),c(Y,Ypred),type='n',xlab='',ylab='',...)
        graphics::points(X,Y)
        graphics::lines(X,Ypred)
        graphics::mtext(side=1,text='t',line=2)
        graphics::mtext(side=2,text=ylab,line=2)
    }
    graphics::par(oldpar)
}

SS_b0g <- function(b0g,spot,groups,dc=NULL,mse=FALSE){
    den <- groups$den
    B0G <- b0g2list(b0g=b0g,groups=groups)
    b0 <- B0G$b0
    g <- B0G$g
    nb <- length(b0)
    nele <- names(groups$num)
    Dp <- betapars(spot=spot,ion=den,dc=dc)
    exp_a0D <- get_exp_a0D(b0g=b0g,spot=spot,groups=groups,dc=dc)
    SS <- (Dp$bkg + exp_a0D - Dp$sig)^2
    for (ele in nele){
        num <- groups$num[[ele]]
        ni <- length(num)
        for (i in 1:ni){
            ion <- num[i]
            Np <- betapars(spot=spot,ion=ion,dc=dc)
            bND <- b0[ion] + g[ele]*Dp$t + dc['g',ion]*(Np$t-Dp$t)
            exp_a0N <- exp_a0D*exp(bND)
            SS <- SS + (Np$bkg + exp_a0N - Np$sig)^2
        }
    }
    out <- sum(SS)
    if (mse){
        np <- length(b0g)  # number of parameters
        ne <- length(b0)   # number of equations
        nm <- length(Np$t) # number of measurements per equation
        out <- out/(ne*nm-np)
    }
    out
}

# analytical solution for (D - bkg) where D is the common denominator
get_exp_a0D <- function(b0g,spot,groups,dc=NULL){
    B0G <- b0g2list(b0g=b0g,groups=groups)
    b0 <- B0G$b0
    g <- B0G$g
    num <- names(b0)
    den <- groups$den
    tD <- hours(spot$time[,den])
    NUM <- 0
    DEN <- 0
    for (ion in num){
        bX <- background(spot,ion)
        X <- spot$signal[,ion]
        gX <- dc['g',ion]
        tX <- hours(spot$time[,ion])
        gXD <- g[element(ion)]
        ebXD <- exp( b0[ion] + gXD*tD + gX*(tX-tD) )
        NUM <- NUM + (X-bX)*ebXD
        DEN <- DEN + ebXD^2
    }
    NUM/DEN
}

# splits a pooled logratio slope and intercept vector into two
b0g2list <- function(b0g,groups){
    dele <- element(groups$den)
    nele <- element(names(groups$num))
    if (dele %in% nele) {
        ng <- length(groups$num) - 1
        g <- c(0,tail(b0g,n=ng))
        names(g)[1] <- dele
    } else {
        ng <- length(groups$num)
        g <- tail(b0g,n=ng)
    }
    nb <- length(b0g) - ng
    b0 <- b0g[1:nb]
    list(b0=b0,g=g)
}

# extract data from a spot for logratios calculation
betapars <- function(spot,ion,dc=NULL){
    out <- alphapars(spot=spot,ion=ion)
    if (is.null(dc)) out$g <- 0
    else out$g <- dc['g',ion]
    out
}

# B = output of common_denominator
groupbypairs <- function(B){
    out <- list()
    out$num <- list()
    for (ion in B$num){
        nele <- element(ion)
        if (nele %in% names(out$num)){
            out$num[[nele]] <- c(out$num[[nele]],ion)
        } else {
            out$num[[nele]] <- ion
        }
    }
    out$den <- B$den
    out
}
