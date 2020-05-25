delta <- function(x,log=TRUE){
    out <- x
    logd <- x$lr - rep(x$ref,length(x$snames))
    nd <- length(logd)
    if (log){
        out$d <- 1000*logd
        J <- 1000*diag(nd)
    } else {
        out$d <- 1000*(exp(logd)-1)
        J <- 1000*exp(logd)*diag(nd)
    }
    out$cov <- J %*% x$cov %*% t(J)
    out$ref <- NULL
    class(out) <- append('delta',class(x))
    out
}

data2table <- function(x){
    if (methods::is(x,'stable')){
        out <- addtab(x)
    } else {
        lists <- x
        lists$snames <- NULL # remove the sample names
        out <- NULL
        for (item in lists){
            out <- cbind(out,addtab(dat=item))
        }
    }
    rownames(out) <- x$snames
    out
}

addtab <- function(dat){
    ni <- length(dat$num)
    ns <- length(dat$lr)/ni
    nc <- 2*ni+ni*(ni-1)/2
    cnames <- rep('',nc)
    DP <- paste0(dat$num,'/',dat$den)
    cnames[2*(1:ni)-1] <- DP
    cnames[2*(1:ni)] <- paste0('s[',DP,']')
    if (nc>2*ni){
        for (i in 1:(ni-1)){
            for (j in (i+1):ni){
                cnames[2*ni+(i-1)*ni+(j-1)] <- paste0('r[',DP[i],',',DP[j],']')
            }
        }
    }
    out <- matrix(0,ns,nc)
    colnames(out) <- cnames
    for (i in 1:ns){
        ii <- (i-1)*ni+(1:ni)
        out[i,2*(1:ni)-1] <- exp(dat$lr[ii])
        J <- exp(dat$lr[ii])*diag(ni)
        E <- J %*% dat$cov[ii,ii] %*% t(J)
        out[i,2*(1:ni)] <- sqrt(diag(E))
        cormat <- cov2cor(E)
        if (nc>2*ni) out[i,(2*ni):nc] <- cormat[upper.tri(cormat)]
    }
    out
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
