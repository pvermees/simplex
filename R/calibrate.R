calibrate <- function(cal,exterr=FALSE){
    if (stable(cal)) out <- calibrate_stable(dat=cal,exterr=exterr)
    else out <- calibrate_geochron(dat=cal,exterr=exterr)
    out
}

calibrate_stable <- function(dat,exterr=FALSE){
    out <- list()
    cal <- dat$stand$fetch(dat)
    snames <- names(dat$x)
    ns <- length(snames)
    nr <- length(dat$num)
    out$snames <- snames
    out$num <- dat$num
    out$den <- dat$den
    out$lr <- rep(0,nr*ns)
    E <- matrix(0,nr*(ns+2),nr*(ns+2))
    J <- matrix(0,nr*ns,nr*(ns+2))
    J[1:(ns*nr),1:(ns*nr)] <- diag(ns*nr)
    for (i in 1:ns){
        sname <- snames[i]
        selected <- (i-1)*nr + (1:nr)
        sp <- spot(dat,sname=sname)
        out$lr[selected] <- sp$lr$b0g[1:nr] + cal$lr - dat$cal$x
        E[selected,selected] <- sp$lr$cov[1:nr,1:nr]
        if (exterr){
            J[selected,ns*nr+(1:nr)] <- diag(nr)
            J[selected,(ns+1)*nr+(1:nr)] <- -diag(nr)
        }
    }
    E[ns*nr+(1:nr),ns*nr+(1:nr)] <- cal$cov
    E[(ns+1)*nr+(1:nr),(ns+1)*nr+(1:nr)] <- dat$cal$cov
    out$cov <- J %*% E %*% t(J)
    class(out) <- "calibrated"
    out
}

calibrate_geochron <- function(cal,exterr=FALSE){
    
}
