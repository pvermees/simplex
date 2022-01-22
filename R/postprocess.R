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
    ref <- cd$calibration$stand$ref
    logd <- cd$calibrated$val - rep(ref$val,length(cd$samples))
    nd <- length(logd)
    if (log){
        del$val <- 1000*logd
        J <- 1000*diag(nd)
    } else {
        del$val <- 1000*(exp(logd)-1)
        J <- 1000*exp(logd)*diag(nd)
    }
    del$cov <- J %*% cd$calibrated$cov %*% t(J)
    nms <- rep(paste0('delta(',names(ref$val),')'),length(cd$samples))
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
#' @rdname data2table
#' @export
data2table.delta <- function(x,cov=FALSE,...){
    del <- x$delta
    cal <- x$calibrated
    data2table_helper(val=del$val,E=del$cov,snames=cal$snames,
                      ratios=cal$ratios,cov=cov)
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
            cormat <- cov2cor(E)
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
        colnames(out) <- c('ratios',rep(ratios,ns))
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
        cormat <- cov2cor(E)
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
    if (identical(method,'U-Pb')){
        iratios <- c('Pb204/Pb206','Pb207/Pb206','Pb206/U238')
        oratios <- c('U238/Pb206','Pb207/Pb206','Pb204/Pb206')
        j <- rbind(c(0,0,-1),
                   c(0,1,0),
                   c(1,0,0))
    } else if (identical(method,'U-Th-Pb')){
        iratios <- c('Pb204/Pb206','Pb207/Pb206','Pb204/Pb208','Pb206/U238','Pb208/Th232')
        oratios <- c('U238/Pb206','Pb207/Pb206','Pb208/Pb206','Th232/U238')
        j <- rbind(c(0,0,0,-1,0),
                   c(0,1,0,0,0),
                   c(1,0,-1,0,0),
                   c(1,0,-1,1,-1))
        format <- 8
    } else if (identical(method,'Th-Pb')){
        iratios <- c('Pb204/Pb208','Pb208/Th232')
        iratios <- c('Th232/Pb208','Pb204/Pb208')
        j <- rbind(c(0,-1),
                   c(-1,0))
    }
    if (all(iratios %in% dat$calibrated$ratios)){
        cal <- dat$calibrated
        snames <- cal$snames
        ns <- length(snames)
        ni <- length(iratios)
        no <- length(oratios)
        J <- matrix(0,no*ns,ni*ns)
        for (i in 1:ns){
            i1 <- (i-1)*no+(1:no)
            i2 <- (i-1)*ni+(1:ni)
            J[i1,i2] <- j
        }
        val <- as.vector(exp(J %*% cal$val))
        E <- diag(val) %*% J %*% cal$cov %*% t(J) %*% diag(val)
        tab <- data2table_helper(val=val,E=E,snames=snames,ratios=oratios)
    } else {
        stop("Input ratios should include: ",iratios)
    }
    if (identical(method,'U-Pb')){
        out <- IsoplotR:::as.UPb(tab,format=5)
    } else if (identical(method,'U-Th-Pb')){
        out <- IsoplotR:::as.UPb(tab,format=8)
    } else if (identical(method,'Th-Pb')){
        out <- IsoplotR:::as.ThPb(tab,format=2)
    }
    out
}
