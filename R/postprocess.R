#' @title calculate delta values
#' @description convert stable isotope logratio data to delta notation
#' @param cd an object of class \code{calibrated}
#' @param log use traditional delta notaion (\code{log=FALSE}) or the
#'     logratio definition (\code{log=TRUE})?
#' @return an object of class \code{delta}
#' @examples
#' data('oxygen',package='simplex')
#' dc <- drift(x=oxygen)
#' lr <- logratios(x=dc)
#' cal <- calibration(lr=lr,stand=standard(preset='NBS28'))
#' cd <- calibrate(cal)
#' del <- delta(cd)
#' tab <- data2table(del)
#' @export
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

#' @title export to table
#' @description convert calibrated SIMS data objects to a flat table
#' @param x an object of class \code{calibrated} or \code{delta}
#' @param ... optional arguments
#' @return a matrix
#' @examples
#' data('Cameca',package='simplex')
#' dc <- drift(x=Cameca)
#' lr <- logratios(x=dc)
#' cal <- calibration(lr=lr,stand=standard(preset='Plesovice'))
#' cd <- calibrate(cal)
#' tab <- data2table(cd)
#' @rdname data2table
#' @export
data2table <- function(x,...){ UseMethod("data2table",x) }
#' @rdname data2table
#' @export
data2table.default <- function(x,...){
    data2table_helper(x=x,option='default',...)
}
#' @rdname data2table
#' @export
data2table.calibrated <- function(x,...){
    data2table_helper(x=x,option='calibrated',...)
}
#' @rdname data2table
#' @export
data2table.delta <- function(x,...){
    data2table_helper(x=x,option='delta',...)
}
data2table_helper <- function(x,option,...){
    p <- data2table_pars(x=x,option=option)
    ni <- length(p$num)
    ns <- length(p$lr)/ni
    nc <- 2*ni+ni*(ni-1)/2
    cnames <- rep('',nc)
    DP <- paste0(p$num,'/',p$den)
    cnames[2*(1:ni)-1] <- DP
    cnames[2*(1:ni)] <- paste0('s[',DP,']')
    if (nc>2*ni){
        for (i in 1:(ni-1)){
            for (j in (i+1):ni){
                cnames[2*ni+(i-1)*(ni-1)+(j-i)] <-
                    paste0('r[',DP[i],',',DP[j],']')
            }
        }
    }
    out <- matrix(0,ns,nc)
    colnames(out) <- cnames
    for (i in 1:ns){
        ii <- (i-1)*ni+(1:ni)
        lr <- as.vector(p$J %*% p$lr[ii])
        covmat <- p$J %*% p$cov[ii,ii] %*% t(p$J)
        if (option=='delta'){
            out[i,2*(1:ni)-1] <- lr
            E <- covmat
        } else {
            out[i,2*(1:ni)-1] <- exp(lr)
            J <- diag(exp(lr),nrow=ni,ncol=ni)
            E <- J %*% covmat %*% t(J)
        }
        out[i,2*(1:ni)] <- sqrt(diag(E))
        cormat <- stats::cov2cor(E)
        if (nc>2*ni) out[i,(2*ni+1):nc] <- cormat[upper.tri(cormat)]
    }
    rownames(out) <- names(x$samples)
    out
}

data2table_pars <- function(x,option){
    type <- datatype(x)
    if (type=='oxygen'){
        num <- c('O17','O18')
        den <- c('O16','O16')
        J <- diag(2)
    } else if (type=='U-Pb'){
        num <- c('U238','Pb207','Pb204')
        den <- c('206Pb','Pb206','Pb206')
        J <- matrix(0,3,3)
        J[1,1] <- -1
        J[2,3] <- 1
        J[3,2] <- 1
    } else if (type=='U-Th-Pb'){
        num <- c('U238','Pb207','Pb204','Pb208','Th232')
        den <- c('206Pb','Pb206','Pb206','Pb206','Pb208')
        J <- matrix(0,5,5)
        J[1,1] <- -1
        J[2,3] <- 1
        J[3,2] <- 1
        J[4,4] <- 1
        J[5,5] <- -1
    } else {
        stop('invalid data type')
    }
    if (option == 'calibrated'){
        lr <- x$calibrated$lr
        covmat <- x$calibrated$cov
    } else if (option == 'delta'){
        lr <- x$delta$d
        covmat <- x$delta$cov
    } else {
        lr <- x$lr
        covmat <- x$cov
    }
    list(num=num,den=den,lr=lr,cov=covmat,J=J)
}

#' @title plot delta values
#' @description visualises stable isotope data as error ellipses in
#'     delta notation
#' @param x an object of class \code{delta}
#' @param ... optional arguments to be passed on to the generic
#'     \code{plot} function.
#' @examples
#' data('oxygen',package='simplex')
#' dc <- drift(x=oxygen)
#' lr <- logratios(x=dc)
#' cal <- calibration(lr=lr,stand=standard(preset='NBS28'))
#' cd <- calibrate(cal)
#' del <- delta(cd)
#' plot(del)
#' @method plot delta
#' @export
plot.delta <- function(x,...){
    del <- x$delta
    np <- length(del$num)-1     # number of plot panels
    nr <- ceiling(sqrt(np)) # number of rows
    nc <- ceiling(np/nr)    # number of columns
    oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    for (i in 1:nr){
        for (j in (i+1):max(nc,nr+1)){
            y <- delta2york(d=x,i=i,j=j)
            xlab <- substitute(delta^{a}*b,
                               list(a=isotope(ion=del$num[i]),
                                    b=element(ion=del$num[i])))
            ylab <- substitute(delta^{a}*b,
                               list(a=isotope(ion=del$num[j]),
                                    b=element(ion=del$num[j])))
            IsoplotR::scatterplot(y,...)
            graphics::mtext(xlab,side=1,line=2)
            graphics::mtext(ylab,side=2,line=2)
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
