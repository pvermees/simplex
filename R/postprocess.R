#' @title calculate delta values
#' @description convert stable isotope logratio data to delta notation
#' @param cd an object of class \code{calibrated}
#' @param log use traditional delta notaion (\code{log=FALSE}) or the
#'     logratio definition (\code{log=TRUE})?
#' @return an object of class \code{delta}
#' @examples
#' data('Cameca_oxygen',package='simplex')
#' dc <- drift(x=Cameca_oxygen)
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
    ref <- do.call(cd$standard$fetchfun,args=list(dat=cd))$ref
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
    nms <- rep(paste0('delta(',del$num,')'),length(cd$samples))
    rownames(del$cov) <- nms
    colnames(del$cov) <- nms
    names(del$d) <- nms
    out$delta <- del
    class(out) <- unique(append('delta',class(cd)))
    out
}

#' @title export to table
#' @description convert calibrated SIMS data objects to a flat table
#' @param x an object of class \code{calibrated} or \code{delta}
#' @param ... optional arguments
#' @return a matrix
#' @examples
#' data('Cameca_UPb',package='simplex')
#' dc <- drift(x=Cameca_UPb)
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
    stop('No default version of data2table.')
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
#' @rdname data2table
#' @export
data2table.logratios <- function(x,log=TRUE,t=NULL,...){
    num <- x$method$num
    den <- x$method$den
    samps <- x$samples
    nr <- length(samps)
    nn <- length(num)
    nrho <- nn*(nn-1)/2
    nc <- 2*nn + nrho
    out <- matrix(NA,nrow=nr,ncol=nc)
    cnames <- rep(NA,nc)
    ii <- 1
    for (i in 1:nn){
        ratio1 <- paste0(num[i],'/',den[i])
        r1 <- ifelse(log,paste0('ln[',ratio1,']'),ratio1)
        cnames[2*i-1] <- r1
        cnames[2*i] <- paste0('s(',r1,')')
        for (j in i:nn){
            if (i!=j){
                ratio2 <- paste0(num[j],'/',den[j])
                r2 <- ifelse(log,paste0('ln[',ratio2,']'),ratio2)
                cnames[2*nn+ii] <- paste0('r(',r1,',',r2,')')
                ii <- ii + 1
            }
        }
    }
    rownames(out) <- names(samps)
    colnames(out) <- cnames
    if (is.null(t)) t <- stats::median(samps[[1]]$time)
    J <- cbind(diag(nn),diag(nn)*hours(t))
    for (i in 1:nr){
        lr <- samps[[i]]$lr
        b0 <- lr$b0g[1:nn]
        g <- lr$b0g[(nn+1):(2*nn)]
        lrt <- b0 + g * hours(t)
        covlrt <- J %*% lr$cov %*% t(J)
        if (log) {
            out[i,2*(1:nn)-1] <- lrt
            out[i,2*(1:nn)] <- sqrt(diag(covlrt))
            cormat <- cov2cor(covlrt)
        } else {
            out[i,2*(1:nn)-1] <- exp(lrt)
            Jexp <- diag(nn)*exp(lrt)
            covmat <- Jexp %*% covlrt %*% t(Jexp)
            out[i,2*(1:nn)] <- sqrt(diag(covmat))
            cormat <- cov2cor(covmat)
            
        }
        if (nc>2*(nn)) out[i,(2*nn+1):(nc)] <- cormat[upper.tri(cormat)]
    }
    out
}
data2table_helper <- function(x,option,...){
    p <- data2table_pars(x=x,option=option)
    nin <- ncol(p$J)
    nout <- nrow(p$J)
    ns <- length(p$lr)/nin
    nc <- 2*nout+nout*(nout-1)/2
    cnames <- rep('',nc)
    if (option=='delta'){
        DP <- paste0('d(',p$num,')')
    } else {
        DP <- paste0(p$num,'/',p$den)
    }
    cnames[2*(1:nout)-1] <- DP
    cnames[2*(1:nout)] <- paste0('s[',DP,']')
    ii <- 2*nout
    for (i in 1:nout){
        for (j in i:nout){
            if (i!=j){
                ii <- ii + 1
                cnames[ii] <- paste0('r[',DP[i],',',DP[j],']')
            }
        }
    }
    out <- matrix(0,ns,nc)
    colnames(out) <- cnames
    for (i in 1:ns){
        ii <- (i-1)*nin+(1:nin)
        lr <- as.vector(p$J %*% p$lr[ii])
        covmat <- p$J %*% p$cov[ii,ii] %*% t(p$J)
        if (option=='delta'){
            out[i,2*(1:nout)-1] <- lr
            E <- covmat
        } else {
            out[i,2*(1:nout)-1] <- exp(lr)
            J <- diag(exp(lr),nrow=nout,ncol=nout)
            E <- J %*% covmat %*% t(J)
        }
        out[i,2*(1:nout)] <- sqrt(diag(E))
        cormat <- stats::cov2cor(E)
        out[i,(2*nout+1):nc] <- cormat[upper.tri(cormat)]
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
    } else if (type=='sulphur'){
        num <- c('S33','S34','S36')
        den <- c('S32','S32','S32')
        J <- diag(3)
    } else if (type=='U-Pb'){
        num <- c('U238','Pb207','Pb204')
        den <- c('Pb206','Pb206','Pb206')
        J <- matrix(0,3,3)
        J[1,1] <- -1
        J[2,3] <- 1
        J[3,2] <- 1
    } else if (type=='U-Th-Pb'){
        num <- c('U238','Pb207','Pb204','Pb208','Th232','Th232','Pb204')
        den <- c('Pb206','Pb206','Pb206','Pb206','U238','Pb208','Pb208')
        J <- matrix(0,7,5)
        J[1,1] <- -1
        J[2,4] <- 1
        J[3,3] <- 1
        J[4,3] <- 1
        J[4,5] <- -1
        J[5,1] <- 1
        J[5,2] <- -1
        J[5,3] <- 1
        J[5,5] <- -1
        J[6,2] <- -1
        J[7,5] <- 1
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
#' data('Cameca_oxygen',package='simplex')
#' dc <- drift(x=Cameca_oxygen)
#' lr <- logratios(x=dc)
#' cal <- calibration(lr=lr,stand=standard(preset='NBS28'))
#' cd <- calibrate(cal)
#' del <- delta(cd)
#' plot(del)
#' @method plot delta
#' @export
plot.delta <- function(x,...){
    del <- x$delta
    nn <- length(del$num)
    np <- nn*(nn-1)/2       # number of plot panels
    nc <- ceiling(sqrt(np)) # number of rows
    nr <- ceiling(np/nc)    # number of columns
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
    graphics::par(oldpar)
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

# TEMPORARY SOLUTION FOR PAPER. NEEDS TIDYING UP AND GENERALING
age <- function(cd){
    l38 <- IsoplotR::settings('lambda','U238')[1]
    num <- cd$calibrated$num
    den <- cd$calibrated$den
    ni <- length(num)
    ns <- length(cd$samples)
    iPbU <- which(num %in% 'Pb206' & den %in% 'U238')
    i <- (1:ns)*ni-ni+iPbU
    out <- list()
    elr <- exp(cd$calibrated$lr[i])
    out$t68 <- log(elr+1)/l38
    J <- diag(elr/(elr+1))/l38
    out$E68 <- J %*% cd$calibrated$cov[i,i] %*% t(J)
    snames <- names(cd$samples)
    names(out$t68) <- snames
    rownames(out$E68) <- snames
    colnames(out$E68) <- snames
    out
}

#' @title convert to IsoplotR
#' @description convert U-Pb or U-Th-Pb data to an IsoplotR object
#' @param dat an object of class \code{calibrated}
#' @param method sets the format of the \code{IsoplotR} oboject. If
#'     \code{'U-Pb'}, produces a format 5 object of class \code{UPb};
#'     if \code{'U-Th-Pb'}, produces a format 8 object of class
#'     \code{UPb}; if \code{'Th-Pb'}, produces a format 2 object of
#'     class \code{ThPb}.
#' @return an object of class \code{UPb} or \code{ThPb}
#' @examples
#' \dontrun{
#' m <- method('GA-UPb')
#' s <- standard(preset="Temora",prefix='TEM')
#' fname <- system.file('SHRIMP.pd',package='simplex')
#' cd <- process(f=fname,m=m,stand=s)
#' UPb <- simplex2IsoplotR(cd)
#' IsoplotR::concordia(UPb)
#' }
#' @export
simplex2IsoplotR <- function(dat,method='U-Pb'){
    tab <- data2table(dat)
    dt <- datatype(dat)
    if (identical(method,'U-Pb')){
        if (identical(dt,'U-Th-Pb')){
            cols <- c(1:6,15:16,21)
            tab <- tab[,cols,drop=FALSE]
        }
        out <- IsoplotR:::as.UPb(tab,format=5)
    } else if (identical(method,'U-Th-Pb')){
        if (!identical(dt,'U-Th-Pb'))
            stop('Invalid data type or U-Th-Pb dating.')
        cols <- c(1:4,7:10,15,17:18,22:23,30)
        tab <- tab[,cols,drop=FALSE]
        out <- IsoplotR:::as.UPb(tab,format=8)
    } else if (identical(method,'Th-Pb')){
        if (!identical(dt,'U-Th-Pb'))
            stop('Invalid data type or U-Th-Pb dating.')
        cols <- c(11:14,33)
        out <- IsoplotR:::as.ThPb(tab,format=2)
    }
    out
}
