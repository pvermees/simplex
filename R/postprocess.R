delta <- function(cd,log=TRUE){
    out <- list()
    ref <- cd$standard$fetch(cd)$ref
    logd <- cd$calibrated$lr - rep(ref,length(cd$samples))
    nd <- length(logd)
    if (log){
        out$d <- 1000*logd
        J <- 1000*diag(nd)
    } else {
        out$d <- 1000*(exp(logd)-1)
        J <- 1000*exp(logd)*diag(nd)
    }
    out$cov <- J %*% cd$calibrated$cov %*% t(J)
    out$ref <- NULL
    class(out) <- 'delta'
    out
}

data2table <- function(x,...){ UseMethod("data2table",x) }
data2table.default <- function(x,...){
    ni <- length(x$num)
    ns <- length(x$lr)/ni
    nc <- 2*ni+ni*(ni-1)/2
    cnames <- rep('',nc)
    DP <- paste0(x$num,'/',x$den)
    cnames[2*(1:ni)-1] <- DP
    cnames[2*(1:ni)] <- paste0('s[',DP,']')
    if (nc>2*ni){
        for (i in 1:(ni-1)){
            for (j in (i+1):ni){
                cnames[2*ni+(i-1)*(ni-1)+(j-i)] <- paste0('r[',DP[i],',',DP[j],']')
            }
        }
    }
    out <- matrix(0,ns,nc)
    colnames(out) <- cnames
    for (i in 1:ns){
        ii <- (i-1)*ni+(1:ni)
        out[i,2*(1:ni)-1] <- exp(x$lr[ii])
        J <- exp(x$lr[ii])*diag(ni)
        E <- J %*% x$cov[ii,ii] %*% t(J)
        out[i,2*(1:ni)] <- sqrt(diag(E))
        cormat <- cov2cor(E)
        if (nc>2*ni) out[i,(2*ni+1):nc] <- cormat[upper.tri(cormat)]
    }
    out
}
data2table.calibrated <- function(x,...){
    data2table.default(x$calibrated)
}

plot.delta <- function(d,...){
    np <- length(d$num)-1     # number of plot panels
    nr <- ceiling(sqrt(np)) # number of rows
    nc <- ceiling(np/nr)    # number of columns
    oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    for (i in 1:nr){
        for (j in (i+1):max(nc,nr+1)){
            y <- delta2york(d=d,i=i,j=j)
            xlab <- substitute(delta^{a}*b,
                               list(a=isotope(ion=d$num[i]),
                                    b=element(ion=d$num[i])))
            ylab <- substitute(delta^{a}*b,
                               list(a=isotope(ion=d$num[j]),
                                    b=element(ion=d$num[j])))
            IsoplotR::scatterplot(y,...)
            mtext(xlab,side=1,line=2)
            mtext(ylab,side=2,line=2)
        }
    }
}

delta2york <- function(d,i,j){
    ns <- length(d$snames)
    ni <- length(d$num)
    ii <- (1:ns)*ni-ni+i
    jj <- (1:ns)*ni-ni+j
    X <- d$d[ii]
    Y <- d$d[jj]
    sX <- sqrt(diag(d$cov)[ii])
    sY <- sqrt(diag(d$cov)[jj])
    rXY <- d$cov[cbind(ii,jj)]/(sX*sY)
    out <- cbind(X,sX,Y,sY,rXY)
    colnames(out) <- c('X','sX','Y','sY','rXY')
    rownames(out) <- d$snames
    out
}
