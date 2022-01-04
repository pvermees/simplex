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
data2table.calibrated <- function(x,snames=NULL,i=NULL,
                                  ratios=NULL,j=NULL,
                                  cov=FALSE,log=TRUE,...){
    cal <- x$calibrated
    if (!is.null(snames)){
        i <- which(cal$snames %in% snames)
    } else if (is.null(i)){
        i <- 1:length(cal$snames)
    }
    if (!is.null(ratios)){
        j <- which(cal$ratios %in% ratios)
    } else if (is.null(j)){
        j <- 1:length(cal$ratios)
    }
    ns <- length(i)
    nr <- length(j)
    ii <- rep((i-1)*nr,each=nr) + rep(j,ns)
    if (log){
        val <- cal$val[ii]
        E <- cal$cov[ii,ii]
    } else {
        val <- exp(cal$val[ii])
        E <- diag(val) %*% cal$cov[ii,ii] %*% diag(val)
    }
    if (cov){
        out <- cbind(val,E)
        colnames(out) <- c('ratios',rep(cal$ratios[j],ns))
        rownames(out) <- rep(cal$snames[i],each=nr)
    } else {
        nc <- 2*nr+nr*(nr-1)/2
        out <- matrix(0,nrow=ns,ncol=nc)
        rownames(out) <- cal$snames[i]
        cnames <- rep(NA,nc)
        out[,2*(1:nr)-1] <- matrix(val,nrow=ns,ncol=nr,byrow=TRUE)
        out[,2*(1:nr)] <- matrix(sqrt(diag(E)),nrow=ns,ncol=nr,byrow=TRUE)
        cnames[2*(1:nr)-1] <- cal$ratios[j]
        cnames[2*(1:nr)] <- paste0('s[',cal$ratios[j],']')
        ci <- 2*nr
        cormat <- cov2cor(E)
        if (nr>1){
            for (r1 in 1:(nr-1)){
                i1 <- (0:(ns-1))*nr + r1
                for (r2 in (r1+1):nr){
                    ci <- ci + 1
                    i2 <- (0:(ns-1))*nr + r2
                    out[,ci] <- diag(cormat[i1,i2])
                    cnames[ci] <- paste0('r[',cal$ratios[j][r1],
                                         ',',cal$ratios[j][r2],']')
                }
            }
        }
        colnames(out) <- cnames
    }
    out
}
#' @rdname data2table
#' @export
data2table.delta <- function(x,...){
    data2table_helper(x=x,option='delta',...)
}
#' @rdname data2table
#' @export
data2table.logratios <- function(x,log=TRUE,t=NULL,addxy=FALSE,...){
    if (addxy & identical(x$method$instrument,'SHRIMP')){
        stop('SHRIMP data do not contain x-y stage coordinates.')
    }
    tavg <- time_average(x,t=t)
    ns <- length(tavg) # number of aliquots
    nr <- length(tavg[[1]]$val) # number of ratios
    si <- ifelse(addxy,2,0) # start index
    nc <- si + 2*nr + nr*(nr-1)/2 # number of columns
    out <- matrix(NA,nrow=ns,nc=nc)
    cnames <- rep(NA,nc)
    ratios <- paste0(x$method$num,'/',x$method$den)
    comb <- utils::combn(ratios,m=2)
    if (log){
        cnames[si+2*(1:nr)-1] <- paste0('ln[',ratios,']')
        cnames[si+2*(1:nr)] <- paste0('s(ln[',ratios,'])')
        if (nr>1) cnames[(si+2*nr+1):nc] <-
                      paste0('r(ln[',comb[1,],'],ln[',comb[2,],'])')
    } else {
        cnames[si+2*(1:nr)-1] <- ratios
        cnames[si+2*(1:nr)] <- paste0('s(',ratios,')')
        if (nr>1) cnames[(si+2*nr+1):nc] <-
                      paste0('r(',comb[1,],',',comb[2,],')')
    }
    if (addxy) cnames[1:2] <- c('x','y')
    colnames(out) <- cnames
    rownames(out) <- names(x$samples)
    for (i in 1:ns){
        if (log){
            val <- tavg[[i]]$val
            E <- tavg[[i]]$cov
        } else {
            val <- exp(tavg[[i]]$val)
            J <- diag(nr) * as.vector(val)
            E <- J %*% tavg[[i]]$cov %*% t(J)
        }
        if (addxy){
            out[i,1] <- x$samples[[i]]$x
            out[i,2] <- x$samples[[i]]$y
        }
        out[i,si+2*(1:nr)-1] <- val
        out[i,si+2*(1:nr)] <- sqrt(diag(E))
        if (nr>1){
            cormat <- cov2cor(E)
            out[i,(si+2*nr+1):nc] <- cormat[upper.tri(cormat)]
        }
    }
    out
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
