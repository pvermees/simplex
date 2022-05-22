#' @title calculate delta values
#' @description convert stable isotope logratio data to delta notation
#' @param cd an object of class \code{calibrated}
#' @param ref the logratio composition of a reference material, i.e. a
#'     list with items \code{val} (vector of logratios) and \code{cov}
#'     (the covariance matrix of \code{val})
#' @param log use traditional delta notation (\code{log=FALSE}) or the
#'     logratio definition (\code{log=TRUE})?
#' @return an object of class \code{delta}
#' @examples
#' \dontrun{
#' data('Cameca_oxygen',package='simplex')
#' dc <- drift(x=Cameca_oxygen)
#' lr <- logratios(x=dc)
#' cal <- calibration(lr=lr,stand=standard(preset='NBS28-O'))
#' cd <- calibrate(cal)
#' del <- delta(cd)
#' tab <- data2table(del)
#' }
#' @export
delta <- function(cd,ref,log=TRUE){
    out <- cd
    del <- list()
    if (missing(ref)){
        ref <- cd$calibration$stand$ref
    }
    ratios <- cd$calibrated$ratios
    logd <- cd$calibrated$val - rep(ref$val[ratios],length(cd$samples))
    nd <- length(logd)
    if (log){
        del$val <- 1000*logd
        J <- 1000*diag(nd)
    } else {
        del$val <- 1000*(exp(logd)-1)
        J <- 1000*exp(logd)*diag(nd)
    }
    del$cov <- J %*% cd$calibrated$cov %*% t(J)
    nms <- rep(paste0('delta(',names(ref$val[ratios]),')'),length(cd$samples))
    rownames(del$cov) <- colnames(del$cov) <- names(del$val) <- nms
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
#' \dontrun{
#' data('Cameca_UPb',package='simplex')
#' dc <- drift(x=Cameca_UPb)
#' lr <- logratios(x=dc)
#' cal <- calibration(lr=lr,stand=standard(preset='Plesovice-t'))
#' cd <- calibrate(cal)
#' tab <- data2table(cd)
#' }
#' @rdname data2table
#' @export
data2table <- function(x,...){ UseMethod("data2table",x) }
#' @rdname data2table
#' @export
data2table.default <- function(x,...){
    stop('No default version of data2table.')
}

#' @param cov logical. If \code{TRUE}, produces an output table of
#'     size \eqn{nm\times{(nm+1)}}, where \eqn{n} is the number of
#'     aliquots and \eqn{m} is the number of logratios. The first
#'     column then represents the concatenated vector of logratios for
#'     all the samples, and the subsequent \eqn{nm} columns contain
#'     its covariance matrix, including inter-sample error
#'     correlations.
#' @param log logical. If \code{TRUE}, returns logratios.
#' @param log4lab logical. If \code{TRUE}, adds `ln[' to the column
#'     labels if \code{log=TRUE}.
#' @rdname data2table
#' @export
data2table.calibrated <- function(x,cov=FALSE,log=TRUE,log4lab=TRUE,...){
    cal <- x$calibrated
    ratios <- cal$ratios
    if (log){
        val <- cal$val
        E <- cal$cov
        if (log4lab) ratios <- paste0('ln(',ratios,')')
    } else {
        val <- exp(cal$val)
        E <- diag(val) %*% cal$cov %*% diag(val)
    }
    data2table_helper(val=val,E=E,snames=cal$snames,
                      ratios=ratios,cov=cov)
}
#' @param cov logical. If \code{TRUE}, returns the logratios along
#'     with the full covariance matrix including inter-sample error
#'     correlations.
#' @rdname data2table
#' @export
data2table.delta <- function(x,cov=FALSE,...){
    del <- x$delta
    cal <- x$calibrated
    data2table_helper(val=del$val,E=del$cov,snames=cal$snames,
                      ratios=cal$ratios,cov=cov)
}

#' @param t The time to which the logratio signal should be
#'     interpolated.
#' @param addxy logical. If \code{TRUE}, adds two columns to the
#'     output table with the x- and y-positions of the different
#'     aliquots. Only relevant to Cameca data, so ignored for SHRIMP
#'     data.
#' @rdname data2table
#' @export
data2table.logratios <- function(x,log=TRUE,t=NULL,addxy=FALSE,...){
    if (addxy & identical(x$method$instrument,'SHRIMP')){
        warning('SHRIMP data do not contain x-y stage coordinates.')
        addxy <- FALSE
    }
    tavg <- time_average(x,t=t)
    ns <- length(tavg) # number of aliquots
    nr <- length(tavg[[1]]$val) # number of ratios
    si <- ifelse(addxy,2,0) # start index
    nc <- si + 2*nr + nr*(nr-1)/2 # number of columns
    out <- matrix(NA,nrow=ns,ncol=nc)
    cnames <- rep(NA,nc)
    ratios <- paste0(x$method$num,'/',x$method$den)
    if (nr>1) comb <- utils::combn(ratios,m=2)
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
            cormat <- stats::cov2cor(E)
            out[i,(si+2*nr+1):nc] <- cormat[upper.tri(cormat)]
        }
    }
    out
}
data2table_helper <- function(val,E,snames,ratios,cov=FALSE){
    ns <- length(snames)
    nr <- length(ratios)
    if (cov){
        out <- cbind(val,E)
        colnames(out) <- c('logratios',rep(ratios,ns))
        rownames(out) <- rep(snames,each=nr)
    } else {
        nc <- 2*nr+nr*(nr-1)/2
        out <- matrix(0,nrow=ns,ncol=nc)
        rownames(out) <- snames
        out[,2*(1:nr)-1] <- matrix(val,nrow=ns,ncol=nr,byrow=TRUE)
        out[,2*(1:nr)] <- matrix(sqrt(diag(E)),nrow=ns,ncol=nr,byrow=TRUE)
        cnames <- rep(NA,nc)
        cnames[2*(1:nr)-1] <- ratios
        cnames[2*(1:nr)] <- paste0('s[',ratios,']')
        ci <- 2*nr
        cormat <- stats::cov2cor(E)
        if (nr>1){
            for (r1 in 1:(nr-1)){
                i1 <- (0:(ns-1))*nr + r1
                for (r2 in (r1+1):nr){
                    ci <- ci + 1
                    i2 <- (0:(ns-1))*nr + r2
                    out[,ci] <- diag(cormat[i1,i2])
                    cnames[ci] <- paste0('r[',ratios[r1],',',ratios[r2],']')
                }
            }
        }
        colnames(out) <- cnames
    }
    out
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
    ratios <- dat$calibrated$ratios
    if (identical(method,'U-Pb')){
        iratios <- c('Pb204/Pb206','Pb207/Pb206','Pb206/U238')
        oratios <- c('U238/Pb206','Pb207/Pb206','Pb204/Pb206')
    } else if (identical(method,'U-Th-Pb')){
        iratios <- c('Pb204/Pb206','Pb207/Pb206','Pb204/Pb208',
                     'Pb206/U238','Pb208/Th232')
        oratios <- c('U238/Pb206','Pb207/Pb206',
                     'Pb208/Pb206','Th232/U238')
    } else if (identical(method,'Th-Pb')){
        iratios <- c('Pb204/Pb208','Pb208/Th232')
        oratios <- c('Th232/Pb208','Pb204/Pb208')
    }
    if (!all(iratios %in% ratios)){
        stop("Input ratios should include: ",iratios)
    }
    nr <- length(ratios)
    no <- length(oratios)
    j <- matrix(0,no,nr)
    colnames(j) <- ratios
    rownames(j) <- oratios
    if (identical(method,'U-Pb')){
        j['U238/Pb206','Pb206/U238'] <- -1
        j['Pb207/Pb206','Pb207/Pb206'] <- 1
        j['Pb204/Pb206','Pb204/Pb206'] <- 1
    } else if (identical(method,'U-Th-Pb')){
        j['U238/Pb206','Pb206/U238'] <- -1
        j['Pb207/Pb206','Pb207/Pb206'] <- 1
        j['Pb208/Pb206','Pb204/Pb206'] <- 1
        j['Pb208/Pb206','Pb204/Pb208'] <- -1
        j['Th232/U238','Pb204/Pb206'] <- 1
        j['Th232/U238','Pb204/Pb208'] <- -1
        j['Th232/U238','Pb206/U238'] <- 1
        j['Th232/U238','Pb208/Th232'] <- -1
    } else if (identical(method,'Th-Pb')){
        j['Th232/Pb208','Pb208/Th232'] <- -1
        j['Pb204/Pb208','Pb204/Pb208'] <- 1
    }
    cal <- dat$calibrated
    snames <- cal$snames
    ns <- length(snames)
    J <- matrix(0,no*ns,nr*ns)
    for (i in 1:ns){
        i1 <- (i-1)*no+(1:no)
        i2 <- (i-1)*nr+(1:nr)
        J[i1,i2] <- j
    }
    val <- as.vector(exp(J %*% cal$val))
    E <- diag(val) %*% J %*% cal$cov %*% t(J) %*% diag(val)
    tab <- data2table_helper(val=val,E=E,snames=snames,ratios=oratios)
    if (identical(method,'U-Pb')){
        out <- IsoplotR:::as.UPb(tab,format=5)
    } else if (identical(method,'U-Th-Pb')){
        out <- IsoplotR:::as.UPb(tab,format=8)
    } else if (identical(method,'Th-Pb')){
        out <- IsoplotR:::as.ThPb(tab,format=2)
    }
    out
}
