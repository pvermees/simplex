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

calibrate_geochron <- function(dat,exterr=FALSE){
    out <- list()
    out$snames <- names(dat$x)
    type <- datatype(dat)
    if (type=="U-Pb"){
        out$PbU <- fractical(dat,type="U-Pb",exterr=exterr)
        out$PbPb <- nofractical(dat,type="U-Pb")
    } else if (type=="U-Th-Pb"){
        out$PbU <- fractical(dat,type="U-Pb",exterr=exterr)
        out$PbPb <- nofractical(dat,type="U-Th-Pb")
        out$ThU <- fractical(dat,type="Th-Pb",exterr=exterr)
    } else {
        stop("Invalid data type supplied to calibrate function.")
    }
    out
}

fractical <- function(dat,type="U-Pb",exterr=FALSE){
    cal <- dat$stand$fetch(dat)
    if (type=='U-Pb'){
        num='Pb206'
        den='U238'
    } else if (type=='Th-Pb'){
        num='Pb208'
        den='Th232'
    }
    snames <- names(dat$x)
    ns <- length(snames)
    out <- list()
    out$num <- num
    out$den <- den
    fit <- dat$cal[[type]]$fit
    fitcov <- diag(c(fit$a[2],fit$b[2]))^2
    fitcov[1,2] <- fit$cov.ab/(fit$a[2]*fit$b[2])
    DP <- paste0(num,den)
    E <- matrix(0,ns+3,ns+3)
    E[ns+1,ns+1] <- cal$cov[DP,DP]
    E[ns+(2:3),ns+(2:3)] <- fitcov
    J <- matrix(0,ns,ns+3)
    out$lr <- rep(0,ns)
    rownames(J) <- snames
    names(out$lr) <- snames
    for (i in 1:ns){
        sp <- spot(dat,i=i)
        b0g <- sp$lr$b0g
        bXlab <- paste0('b0[',sp$cal[[type]]$oxide,'/',den,']')
        gXlab <- paste0('g[',sp$cal[[type]]$oxide,'/',den,']')
        bYlab <- paste0('b0[',num,'/',den,']')
        gYlab <- paste0('g[',num,'/',den,']')
        tt <- sp$cal$t
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
        out$lr[i] <- XY[2] + cal$lr[DP] - (fit$a[1]+fit$b[1]*XY[1])
        J[i,i] <- -fit$b[1]  # dlrdX
        J[i,i+1] <- 1        # dlrdY
        if (exterr){
            J[i,ns+1] <- 1       # dlrdcal68
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
    snames <- names(dat$x)
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
