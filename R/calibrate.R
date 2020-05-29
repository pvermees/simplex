#' @title calibrate SIMS data
#' @description convert signal logratios to atomic logratios
#' @param cal an object of class \code{calibration}
#' @param exterr logical flag indicating whether the uncertainty of
#'     the calibration should be propagated.
#' @return an object of class \code{calibrated}
#' @examples
#' data('oxygen')
#' dc <- drift(x=oxygen)
#' lr <- logratios(x=dc)
#' cal <- calibration(lr=lr,stand=standard(preset='NBS28'))
#' cd <- calibrate(cal)
#' del <- delta(cd,log=FALSE)
#' tab <- data2table(del)
#' @export
calibrate <- function(cal,exterr=FALSE){
    if (stable(cal)) out <- calibrate_stable(dat=cal,exterr=exterr)
    else out <- calibrate_geochron(dat=cal,exterr=exterr)
    out
}

calibrate_stable <- function(dat,exterr=FALSE){
    out <- dat
    scal <- dat$standard$fetch(dat)
    dcal <- dat$calibration
    num <- dat$m$num
    den <- dat$m$den
    snames <- names(dat$samples)
    ns <- length(snames)
    nr <- length(num)
    calibrated <- list()
    calibrated$num <- num
    calibrated$den <- den
    calibrated$lr <- rep(0,nr*ns)
    E <- matrix(0,nr*(ns+2),nr*(ns+2))
    J <- matrix(0,nr*ns,nr*(ns+2))
    J[1:(ns*nr),1:(ns*nr)] <- diag(ns*nr)
    for (i in 1:ns){
        sname <- snames[i]
        selected <- (i-1)*nr + (1:nr)
        sp <- spot(dat,sname=sname)
        calibrated$lr[selected] <- sp$lr$b0g[1:nr] + scal$lr - dcal$lr
        E[selected,selected] <- sp$lr$cov[1:nr,1:nr]
        if (exterr){
            J[selected,ns*nr+(1:nr)] <- diag(nr)
            J[selected,(ns+1)*nr+(1:nr)] <- -diag(nr)
        }
    }
    E[ns*nr+(1:nr),ns*nr+(1:nr)] <- scal$cov
    E[(ns+1)*nr+(1:nr),(ns+1)*nr+(1:nr)] <- dcal$cov
    calibrated$cov <- J %*% E %*% t(J)
    out$calibrated <- calibrated
    class(out) <- append("calibrated",class(out))
    out
}

calibrate_geochron <- function(dat,exterr=FALSE){
    out <- dat
    type <- datatype(dat)
    if (type=="U-Pb"){
        UPb <- fractical(dat,type="U-Pb",exterr=exterr)
        PbPb <- nofractical(dat,type="U-Pb")
        out$calibrated <- mergecal(UPb,PbPb)
    } else if (type=="U-Th-Pb"){
        UPb <- fractical(dat,type="U-Pb",exterr=exterr)
        PbPb <- nofractical(dat,type="U-Th-Pb")
        ThPb <- fractical(dat,type="Th-Pb",exterr=exterr)
        out$calibrated <- mergecal(UPb,ThPb,PbPb)
    } else {
        stop("Invalid data type supplied to calibrate function.")
    }
    class(out) <- append("calibrated",class(out))
    out
}

fractical <- function(dat,type="U-Pb",exterr=FALSE){
    scal <- dat$standard$fetch(dat) # standard calibration
    dcal <- dat$calibration[[type]]      # data calibration
    if (type=='U-Pb'){
        num='Pb206'
        den='U238'
    } else if (type=='Th-Pb'){
        num='Pb208'
        den='Th232'
    }
    snames <- names(dat$samples)
    ns <- length(snames)
    out <- list()
    out$num <- num
    out$den <- den
    fit <- dcal$fit
    fitcov <- matrix(c(fit$a[2]^2,fit$cov.ab,
                       fit$cov.ab,fit$b[2]^2),2,2)
    DP <- paste0(num,den)
    E <- matrix(0,ns+3,ns+3)
    E[ns+1,ns+1] <- scal$cov[DP,DP]
    E[ns+(2:3),ns+(2:3)] <- fitcov
    J <- matrix(0,ns,ns+3)
    out$lr <- rep(0,ns)
    rownames(J) <- snames
    names(out$lr) <- snames
    for (i in 1:ns){
        sp <- spot(dat,i=i)
        b0g <- sp$lr$b0g
        bXlab <- paste0('b0[',dcal$oxide,'/',den,']')
        gXlab <- paste0('g[',dcal$oxide,'/',den,']')
        bYlab <- paste0('b0[',num,'/',den,']')
        gYlab <- paste0('g[',num,'/',den,']')
        tt <- dcal$t
        XY <- rep(0,2)
        XY[1] <- b0g[bXlab] + tt*b0g[gXlab]
        XY[2] <- b0g[bYlab] + tt*b0g[gYlab]
        Eb0g <- sp$lr$cov[c(bXlab,gXlab,bYlab,gYlab),
                          c(bXlab,gXlab,bYlab,gYlab)]
        Jb0g <- matrix(0,2,4)
        Jb0g[1,1] <- 1
        Jb0g[1,2] <- tt
        Jb0g[2,3] <- 1
        Jb0g[2,4] <- tt
        E[i+(0:1),i+(0:1)] <- Jb0g %*% Eb0g %*% t(Jb0g)
        out$lr[i] <- XY[2] + scal$lr[DP] - (fit$a[1]+fit$b[1]*XY[1])
        J[i,i] <- -fit$b[1]  # dlrdX
        J[i,i+1] <- 1        # dlrdY
        if (exterr){
            J[i,ns+1] <- 1       # dlrdscalDP
            J[i,ns+2] <- -1      # dlrda
            J[i,ns+3] <- -XY[1]  # dlrdb
        }    
    }
    out$cov <- J %*% E %*% t(J)
    out
}

nofractical <- function(dat,type="U-Pb"){
    if (type=='U-Pb'){
        num=c('Pb204','Pb207')
        den=c('Pb206','Pb206')
    } else if (type=='U-Th-Pb'){
        num=c('Pb204','Pb207','Pb208')
        den=c('Pb206','Pb206','Pb206')
    }    
    snames <- names(dat$samples)
    ns <- length(snames)
    nr <- length(num)
    out <- list()
    out$num <- num
    out$den <- den
    out$lr <- rep(0,nr*ns)
    out$cov <- matrix(0,nr*ns,nr*ns)
    for (i in 1:ns){
        sp <- spot(dat,i=i)
        j <- (i-1)*nr
        lr <- paste0('b0[',num,'/',den,']')
        out$lr[j+(1:nr)] <- sp$lr$b0g[lr]
        out$cov[j+(1:nr),j+(1:nr)] <- sp$lr$cov[lr,lr]
    }
    out
}

mergecal <- function(...){
    cals <- list(...)
    num <- NULL
    den <- NULL
    for (cal in cals){
        num <- c(num,cal$num)
        den <- c(den,cal$den)
    }
    ni <- length(num)
    ns <- length(cal$lr)/length(cal$num)
    out <- list()
    out$num <- num
    out$den <- den
    out$lr <- rep(0,ni*ns)
    out$cov <- matrix(0,ni*ns,ni*ns)
    for (cal in cals){
        i <- which(num %in% cal$num)
        ii <- as.vector(sapply((0:(ns-1))*ni,'+',i))
        out$lr[ii] <- cal$lr
        out$cov[ii,ii] <- cal$cov
    }
    out
}

#' @title plot calibrated data
#' @description shows the calibrated data on a logratio plot.
#' @param x an object of class \code{calibrated}
#' @param type for U-Pb or U-Th-Pb data, either \code{'U-Pb'} or
#'     \code{'Th-Pb'}.
#' @param ... optional arguments to be passed on to the generic
#'     \code{plot} function
#' @examples
#' data('Cameca')
#' dc <- drift(x=Cameca)
#' lr <- logratios(x=dc)
#' cal <- calibration(lr=lr,stand=standard(preset='Plesovice'))
#' cd <- calibrate(cal)
#' plot(cd,type='U-Pb')
#' @export
plot.calibrated <- function(x,type,...){
    if (missing(type)) cal <- x$calibration[[1]]
    else cal <- x$calibration[[type]]
    y <- beta2york(lr=x,t=cal$t,num=cal$num,den=cal$den)
    IsoplotR::scatterplot(y,...)
}
