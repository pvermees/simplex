delta <- function(cd,log=TRUE){
    out <- cd
    del <- list()
    del$num <- cd$calibrated$num
    del$den <- cd$calibrated$den
    ref <- cd$standard$fetch(cd)$ref
    logd <- cd$calibrated$lr - rep(ref,length(cd$samples))
    nd <- length(logd)
    if (log){
        del$d <- 1000*logd
        J <- 1000*diag(nd)
    } else {
        del$d <- 1000*(exp(logd)-1)
        J <- 1000*exp(logd)*diag(nd)
    }
    del$cov <- J %*% cd$calibrated$cov %*% t(J)
    out$delta <- del
    class(out) <- append('delta',class(cd))
    out
}

data2table <- function(x,...){ UseMethod("data2table",x) }
data2table.default <- function(x,...){
    data2table_helper(x=x,option='default',...)
}
data2table.calibrated <- function(x,...){
    data2table_helper(x=x,option='calibrated',...)
}
data2table.delta <- function(x,...){
    data2table_helper(x=x,option='delta',...)
}
data2table_helper <- function(x,option,...){
    if (option == 'calibrated'){
        X <- x[[option]]
        lr <- 'lr'
    } else if (option == 'delta'){
        X <- x[[option]]
        lr <- 'd'
    } else {
        X <- x
        lr <- 'lr'
    }
    ni <- length(X$num)
    ns <- length(X[[lr]])/ni
    nc <- 2*ni+ni*(ni-1)/2
    cnames <- rep('',nc)
    DP <- paste0(X$num,'/',X$den)
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
        if (option=='delta'){
            out[i,2*(1:ni)-1] <- X[[lr]][ii]
            J <- diag(ni)
        } else {
            out[i,2*(1:ni)-1] <- exp(X[[lr]][ii])
            J <- exp(X[[lr]][ii])*diag(ni)
        }
        E <- J %*% X$cov[ii,ii] %*% t(J)
        out[i,2*(1:ni)] <- sqrt(diag(E))
        cormat <- cov2cor(E)
        if (nc>2*ni) out[i,(2*ni+1):nc] <- cormat[upper.tri(cormat)]
    }
    rownames(out) <- names(x$samples)
    out
}

plot.delta <- function(d,...){
    del <- d$delta
    np <- length(del$num)-1     # number of plot panels
    nr <- ceiling(sqrt(np)) # number of rows
    nc <- ceiling(np/nr)    # number of columns
    oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    for (i in 1:nr){
        for (j in (i+1):max(nc,nr+1)){
            y <- delta2york(d=d,i=i,j=j)
            xlab <- substitute(delta^{a}*b,
                               list(a=isotope(ion=del$num[i]),
                                    b=element(ion=del$num[i])))
            ylab <- substitute(delta^{a}*b,
                               list(a=isotope(ion=del$num[j]),
                                    b=element(ion=del$num[j])))
            IsoplotR::scatterplot(y,...)
            mtext(xlab,side=1,line=2)
            mtext(ylab,side=2,line=2)
        }
    }
}

delta2york <- function(d,i,j){
    del <- d$delta
    ni <- length(del$num)
    ns <- length(del$d)/ni
    ii <- (1:ns)*ni-ni+i
    jj <- (1:ns)*ni-ni+j
    X <- del$d[ii]
    Y <- del$d[jj]
    sX <- sqrt(diag(del$cov)[ii])
    sY <- sqrt(diag(del$cov)[jj])
    rXY <- del$cov[cbind(ii,jj)]/(sX*sY)
    out <- cbind(X,sX,Y,sY,rXY)
    colnames(out) <- c('X','sX','Y','sY','rXY')
    rownames(out) <- names(d$samples)
    out
}
