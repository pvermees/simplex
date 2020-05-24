logratios <- function(x){
    out <- x
    snames <- names(x$x)
    for (sname in snames){
        print(sname)
        sp <- spot(dat=x,sname=sname)
        out$x[[sname]]$lr <- logratios.spot(x=sp)
    }
    class(out) <- append("logratios",class(out))
    out
}

logratios.spot <- function(x){
    num <- x$m$num
    den <- x$m$den
    B <- common_denominator(c(num,den))
    groups <- groupbypairs(B)
    init <- init_logratios(spot=x,groups=groups)
    if (faraday(x)) fn <- faraday_misfit_b0g
    else fn <- sem_misfit_b0g
    fit <- optim(par=init,f=fn,method='L-BFGS-B',
                 lower=init-2,upper=init+2,spot=x,
                 groups=groups,hessian=TRUE)
    fit$cov <- MASS::ginv(fit$hessian)
    pred <- do.call(what=fn,
                    args=list(b0g=fit$par,spot=x,groups=groups,predict=TRUE))
    out <- common2original(fit=fit,num=num,den=den,groups=groups)
    out <- c(out,pred)
    invisible(out)
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
    out <- list()
    out$b0g <- b0gout
    out$cov <- J %*% fit$cov %*% t(J)
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

init_logratios <- function(spot,groups){
    b0 <- NULL
    g <- NULL
    b0names <- NULL
    gnames <- NULL
    den <- groups$den
    for (nums in groups$num){
        if (is.null(spot$dc)){
            for (num in nums){
                Np <- alphapars(spot,num)
                Dp <- alphapars(spot,den)
                absND <- abs((Np$sig-Np$bkg)/(Dp$sig-Dp$bkg))
                b0 <- c(b0,mean(log(absND[absND>0])))
            }
        } else {
            b0 <- c(b0,log(spot$dc['exp_a0',nums]/spot$dc['exp_a0',den]))
        }
        b0names <- c(b0names,nums)
        nele <- unique(element(nums))
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

faraday_misfit_b0g <- function(b0g,spot,groups,predict=FALSE){
    B0G <- b0g2list(b0g=b0g,groups=groups)
    b0 <- B0G$b0
    g <- B0G$g
    D <- betapars(spot=spot,ion=groups$den)
    nele <- length(groups$num)
    obsb <- NULL
    predb <- NULL
    ions <- groups$den
    meas <- (D$sig-D$bkg)
    tt <- D$t
    for (ele in nele){
        num <- groups$num[[ele]]
        ni <- length(num)
        for (i in 1:ni){
            ion <- num[i]
            N <- betapars(spot=spot,ion=ion)
            bND <- b0[ion] + g[ele]*D$t + N$g*(N$t-D$t)
            predb <- cbind(predb,bND)
            obsb <- cbind(obsb, log(N$sig-N$bkg) - log(D$sig-D$bkg))
            ions <- c(ions,ion)
            meas <- cbind(meas,N$sig-N$bkg)
            tt <- cbind(tt,N$t)
        }
    }
    if (predict){
        colnames(meas) <- ions
        colnames(tt) <- ions
        expb <- cbind(1,exp(predb))
        frac <- sweep(expb,1,rowSums(expb),'/')
        colnames(frac) <- ions
        out <- list()
        out$t <- tt
        out$obs <- meas
        out$pred <- sweep(frac,1,rowSums(meas),'*')
        out$outliers <- rep(FALSE,length(D$t))
    } else {
        misfit <- obsb-predb
        covmat <- cov(misfit)
        out <- sum(mahalanobis(x=misfit,center=0*b0,cov=covmat))/2
    }
    out
}

sem_misfit_b0g <- function(b0g,spot,groups,predict=FALSE){
    B0G <- b0g2list(b0g=b0g,groups=groups)
    b0 <- B0G$b0
    g <- B0G$g
    D <- betapars(spot=spot,ion=groups$den)
    pbc <- NULL # predicted count logratios
    ions <- groups$den
    counts <- D$counts
    tt <- D$t
    for (ele in names(groups$num)){
        num <- groups$num[[ele]]
        ni <- length(num)
        for (i in 1:ni){
            ion <- num[i]
            N <- betapars(spot=spot,ion=ion)
            bc <- b0[ion] + g[ele]*D$t + N$g*(N$t-D$t) + log(N$edt) - log(D$edt)
            pbc <- cbind(pbc,bc)
            ions <- c(ions,ion)
            counts <- cbind(counts,N$counts)
            tt <- cbind(tt,N$t)
        }
    }
    thetabkg <- D$bkgcounts/(D$bkgcounts+rowSums(counts))
    nb <- length(b0)
    expbc <- cbind(1-nb*thetabkg,exp(pbc))
    theta <- sweep(sweep(expbc,1,rowSums(expbc),'/'),1,thetabkg)
    colnames(tt) <- ions
    colnames(counts) <- ions
    colnames(theta) <- ions
    if (predict){
        out <- list()
        out$t <- tt
        out$obs <- counts
        out$pred <- sweep(theta,1,rowSums(counts),'*')
        out$outliers <- rep(FALSE,length(D$t))
    } else {
        out <- -sum(log(theta)*counts)
    }
    out
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
betapars <- function(spot,ion){
    out <- alphapars(spot=spot,ion=ion)
    if (!(ion%in%colnames(spot$dc))) out$g <- 0
    else out$g <- spot$dc['g',ion]
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

plot.logratios <- function(x,sname,i=1,option=1,...){
    spot <- spot(x,sname,i=1)
    if (option==1){
        plot_logratios(spot=spot,...)
    } else if (option==2){
        plot_signals(spot=spot,...)
    } else {
        stop("Invalid plot option.")
    }
}

plot_logratios <- function(spot,...){
    bad <- logratios.spot(x=spot)$outliers
    num <- spot$m$num
    den <- spot$m$den
    b0g <- spot$lr$b0g
    nb <- length(b0g)/2
    np <- length(num)       # number of plot panels
    nr <- ceiling(sqrt(np)) # number of rows
    nc <- ceiling(np/nr)    # number of columns
    nt <- nrow(spot$time)
    oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    for (i in 1:np){
        ratio <- paste0(num[i],'/',den[i])
        Np <- betapars(spot=spot,ion=num[i])
        Dp <- betapars(spot=spot,ion=den[i])
        X <- Dp$t
        Y <- (Np$sig-Np$bkg)/(Dp$sig-Dp$bkg)
        b0 <- b0g[paste0('b0[',ratio,']')]
        g <- b0g[paste0('g[',ratio,']')]
        Ypred <- exp(b0 + g*Dp$t + Np$g*(Np$t-Dp$t))
        ylab <- paste0('(',num[i],'-b)/(',den[i],'-b)')
        graphics::plot(c(X,X),c(Y,Ypred),type='n',xlab='',ylab='',...)
        graphics::points(X[!bad],Y[!bad],pch=21)
        graphics::points(X[bad],Y[bad],pch=4)
        graphics::lines(X,Ypred)
        graphics::mtext(side=1,text='t',line=2)
        graphics::mtext(side=2,text=ylab,line=2)
    }
    graphics::par(oldpar)
}

plot_signals <- function(spot,...){
    ions <- spot$m$ions
    np <- length(ions)      # number of plot panels
    nr <- ceiling(sqrt(np)) # number of rows
    nc <- ceiling(np/nr)    # number of columns
    lr <- logratios.spot(x=spot)
    bad <- lr$outliers
    oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    for (ion in ions){
        ylab <- paste0(ion,'- b')
        graphics::plot(lr$t[!bad,ion],lr$obs[!bad,ion],type='n',xlab='',ylab='',...)
        graphics::points(lr$t[!bad,ion],lr$obs[!bad,ion],pch=21)
        graphics::points(lr$t[bad,ion],lr$obs[bad,ion],pch=4)
        graphics::mtext(side=1,text='t',line=2)
        graphics::mtext(side=2,text=ylab,line=2)
        graphics::lines(lr$t[!bad,ion],lr$pred[!bad,ion])
    }
    graphics::par(oldpar)
}
