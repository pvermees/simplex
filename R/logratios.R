#' @title logratio calculation
#' @description extracts logratio averages from time resolved SIMS data
#' @param x an object of class \code{drift}
#' @return an object of class \code{logratios}
#' @examples
#' \dontrun{
#' data('SHRIMP',package='simplex')
#' dc <- drift(x=SHRIMP)
#' lr <- logratios(dc)
#' plot(lr,i=1,option=2)
#' }
#' @export
logratios <- function(x,i=NULL){
    logratios_helper(x,i=i)
}
logratios_helper <- function(x,i=NULL,gui=FALSE){
    out <- x
    snames <- names(x$samples)
    ns <- length(snames)
    if (is.null(i)) ii <- 1:ns
    else ii <- i
    for (j in ii){
        sname <- snames[j]
        if (gui){
            shinylight::sendInfoText(paste(" (processing",sname,")"))
            shinylight::sendProgress(j,ns)
        } else {
            print(sname)
        }
        sp <- spot(dat=x,sname=sname)
        out$samples[[sname]]$lr <- logratios.spot(x=sp)
    }
    class(out) <- unique(append("logratios",class(out)))
    out
}

logratios.spot <- function(x){
    num <- x$method$num
    den <- x$method$den
    B <- common_denominator(c(num,den))
    groups <- groupbypairs(B)
    init <- init_logratios(spot=x,groups=groups)
    if (faraday(x)) fn <- faraday_misfit_b0g
    else fn <- sem_misfit_b0g
    fit <- stats::optim(par=init,f=fn,method='L-BFGS-B',
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
    B0G <- b0g2list(b0g=fit$par,groups=groups,cov=fit$cov)
    b0in <- B0G$b0
    gin <- B0G$g
    bnamesin <- names(b0in)
    gnamesin <- names(gin)
    rnames <- paste0(num,'/',den)
    outnames <- c(paste0('b0[',rnames,']'),
                  paste0('g[',rnames,']'))
    nin <- length(b0in)
    nout <- length(rnames)
    J <- matrix(0,nrow=2*nout,ncol=ncol(B0G$cov))
    colnames(J) <- colnames(B0G$cov)
    rownames(J) <- outnames
    b0gout <- rep(0,2*nout)
    names(b0gout) <- outnames
    for (i in 1:nout){
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
            b0gout[nout+i] <- 0
        } else {
            inele <- which(gnamesin %in% nele)
            b0gout[nout+i] <- b0gout[nin+i] + gin[inele]
            J[nout+i,nin+inele] <- 1
            if (dele != element(groups$den)){
                idele <- which(gnamesin %in% dele)
                b0gout[nout+i] <- b0gout[nout+i] - gin[idele]
                J[nout+i,nin+idele] <- -1
            }
        }
    }
    out <- list()
    out$b0g <- b0gout
    out$cov <- J %*% B0G$cov %*% t(J)
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
            b0 <- c(b0,spot$dc['a0',nums]-spot$dc['a0',den])
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
        for (j in 1:ni){
            ion <- num[j]
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
    } else {
        i <- 1:length(D$t)
        if (!is.null(spot$outliers)) i <- i[-spot$outliers]
        misfit <- obsb[i,,drop=FALSE]-predb[i,,drop=FALSE]
        covmat <- stats::cov(misfit)
        out <- sum(stats::mahalanobis(x=misfit,center=0*b0,cov=covmat))/2
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
    theta <- sweep(sweep(expbc,1,rowSums(expbc),'/'),1,thetabkg,'+')
    colnames(tt) <- ions
    colnames(counts) <- ions
    colnames(theta) <- ions
    if (predict){
        out <- list()
        out$t <- tt
        out$obs <- counts
        out$pred <- sweep(theta,1,rowSums(counts),'*')
    } else {
        i <- 1:length(D$t)
        if (!is.null(spot$outliers)) i <- i[-spot$outliers]
        out <- -stats::dmultinom(counts[i,],prob=theta[i,],log=TRUE)
    }
    out
}

# splits a pooled logratio slope and intercept vector into two
b0g2list <- function(b0g,groups,cov){
    dele <- element(groups$den)
    nele <- element(names(groups$num))
    if (dele %in% nele) {
        ng <- length(groups$num) - 1
        g <- c(0,utils::tail(b0g,n=ng))
        names(g)[1] <- dele
    } else {
        ng <- length(groups$num)
        g <- utils::tail(b0g,n=ng)
    }
    nb <- length(b0g) - ng
    b0 <- b0g[1:nb]
    out <- list(b0=b0,g=g)
    if (!missing(cov)){
        nc <- length(b0) + length(g)
        out$cov <- matrix(0,nc,nc)
        if (ng>0) i <- c(1:nb,nc-ng+(1:ng))
        else i <- 1:nb
        out$cov[i,i] <- cov
        outnames <- c(names(b0),names(g))
        colnames(out$cov) <- outnames
        rownames(out$cov) <- outnames
    }
    out
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

#' @title plot logratios
#' @description plot time resolved logratio data
#' @param x an object of class \code{logratios}
#' @param sname the sample name to be shown
#' @param i the sample number to be shown
#' @param ratios logical. If \code{FALSE}, plots the raw signals
#' versus time. If \code{TRUE}, plots the ratios against time.
#' Both plots show the fitted values as a solid line. Note that, for
#' single collector datasets, the numerator and denominator of the
#' measured ratios correspond to different times.
#' @param ... optional arguments to be passed on to the generic
#'     \code{plot} function.
#' @examples
#' \dontrun{
#' data('SHRIMP',package='simplex')
#' dc <- drift(x=SHRIMP)
#' lr <- logratios(dc)
#' plot(lr,i=1,logratios=TRUE)
#' }
#' @method plot logratios
#' @export
plot.logratios <- function(x,sname=NULL,i=1,ratios=FALSE,...){
    spot <- spot(x,sname=sname,i=i)
    if (ratios){
        plot_ratios(spot=spot,...)
    } else {
        plot_signals(spot=spot,...)
    }
}

plot_ratios <- function(spot,xist=FALSE,...){
    num <- spot$method$num
    den <- spot$method$den
    b0g <- spot$lr$b0g
    nb <- length(b0g)/2
    np <- length(num)       # number of plot panels
    nr <- ceiling(sqrt(np)) # number of rows
    nc <- ceiling(np/nr)    # number of columns
    nt <- nrow(spot$time)
    bg <- rep('black',nt)
    if (!is.null(spot$outliers)) bg[spot$outliers] <- 'white'
    oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    for (i in 1:np){
        ratio <- paste0(num[i],'/',den[i])
        Np <- betapars(spot=spot,ion=num[i])
        Dp <- betapars(spot=spot,ion=den[i])
        if (xist){
            X <- seconds(Dp$t)
            Xlab <- 't (s)'
        } else {
            X <- 1:nt
            Xlab <- 'cycle'
        }
        Y <- (Np$sig-Np$bkg)/(Dp$sig-Dp$bkg)
        b0 <- b0g[paste0('b0[',ratio,']')]
        g <- b0g[paste0('g[',ratio,']')]
        Ypred <- exp(b0 + g*Dp$t + Np$g*(Np$t-Dp$t))
        ylab <- paste0('(',num[i],'-b)/(',den[i],'-b)')
        graphics::plot(c(X,X),c(Y,Ypred),type='n',xlab='',ylab='',...)
        graphics::points(X,Y,pch=21,bg=bg)
        graphics::lines(X,Ypred)
        graphics::mtext(side=1,text=Xlab,line=2)
        graphics::mtext(side=2,text=ylab,line=2)
    }
    graphics::par(oldpar)
}

plot_signals <- function(spot,xist=FALSE,...){
    ions <- colnames(spot$lr$obs)
    np <- length(ions)      # number of plot panels
    nr <- ceiling(sqrt(np)) # number of rows
    nc <- ceiling(np/nr)    # number of columns
    lr <- logratios.spot(x=spot)
    nt <- nrow(spot$time)
    bg <- rep('black',nt)
    if (!is.null(spot$outliers)) bg[spot$outliers] <- 'white'
    oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    for (ion in ions){
        if (xist){
            tt <- seconds(lr$t[,ion])
            tlab <- 't (s)'
        } else {
            tt <- 1:nt
            tlab <- 'cycle'
        }
        graphics::plot(c(tt,tt),c(lr$obs[,ion],lr$pred[,ion]),
                       type='n',xlab='',ylab='',...)
        graphics::points(tt,lr$obs[,ion],pch=21,bg=bg)
        graphics::mtext(side=1,text=tlab,line=2)
        graphics::mtext(side=2,text=ion,line=2)
        graphics::lines(tt,lr$pred[,ion])
    }
    graphics::par(oldpar)
}
