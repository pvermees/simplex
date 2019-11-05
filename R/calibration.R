lr2XY <- function(lr,dat,c64=NULL){
    snames <- names(dat)
    ns <- length(snames)
    out <- list()
    labels <- c('X','Y','L7a','L4a')
    calions <- c('L4c','L6c','L7c','LUc')
    J <- matrix(0,4,4)
    rownames(J) <- labels
    colnames(J) <- calions
    J['X','LUc'] <- 1
    J['Y','L6c'] <- 1
    for (sname in snames){
        d4 <- dat[[sname]]$dwelltime['204Pb']
        d6 <- dat[[sname]]$dwelltime['206Pb']
        d7 <- dat[[sname]]$dwelltime['207Pb']
        dU <- dat[[sname]]$dwelltime['238U']
        dUO <- dat[[sname]]$dwelltime['238U16O2']
        L4c <- lr[[sname]]$x['L4c']
        L6c <- lr[[sname]]$x['L6c']
        L7c <- lr[[sname]]$x['L7c']
        LUc <- lr[[sname]]$x['LUc']
        X <- LUc + log(dU/dUO)
        if (is.null(c64)){
            Y <- L6c + log(dU/d6)
        } else {
            Y <- log(dU) + log(dU/d6) + log(1-(c64*d6)/(d4*exp(L4c)))
            J['Y','L4c'] <- c64*d6/(c64*d6 + exp(L4c)*d4)
        }
        L7a <- L7c + log(d6/d7)
        L4a <- L4c + log(d6/d4)
        J['L7a','L7c'] <- log(d6/d7)
        J['L4a','L4c'] <- log(d6/d4)
        out[[sname]] <- list()
        out[[sname]]$x <- c(X,Y,L7a,L4a)
        names(out[[sname]]$x) <- labels
        out[[sname]]$cov <- J %*% lr[[sname]]$cov[calions,calions] %*% t(J)
        out[[sname]]$omega <- solve(out[[sname]]$cov)
    }
    out
}

# build a calibration curve
calibration <- function(lr,dat,plot=TRUE,disp=TRUE,omit=NULL,...){
    XY <- lr2XY(lr=lr,dat=dat)
    fit <- yorkfit(XY,omit=omit)
    out <- list()
    out$AB <- fit$par
    out$df <- length(dat)-2
    out$mswd <- 2*fit$value/out$df
    if (disp) d <- out$mswd
    else d <- 1
    out$cov <- d*solve(fit$hessian)
    labels <- c('A','B')
    names(out$AB) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    if (plot) calplot(XY,fit=out,omit=omit)
    out
}

yorkfit <- function(XY,omit=NULL){
    if (!is.null(omit)){
        keep <- !((1:length(XY))%in%omit)
        XY <- subset(XY,subset=keep)
    }
    snames <- names(XY)
    ns <- length(snames)
    X <- rep(0,ns)
    Y <- rep(0,ns)
    for (i in 1:ns){
        X[i] <- XY[[i]]$x['X']
        Y[i] <- XY[[i]]$x['Y']
    }
    init <- stats::lm(Y ~ X)$coef
    stats::optim(init,misfit_york,method='BFGS',XY=XY,hessian=TRUE)
}

misfit_york <- function(AB,XY=XY){
    A <- AB[1]
    B <- AB[2]
    SS <- 0
    snames <- names(XY)
    D <- matrix(0,1,2)
    for (sname in snames){
        X <- XY[[sname]]$x['X']
        Y <- XY[[sname]]$x['Y']
        O <- XY[[sname]]$omega[c('X','Y'),c('X','Y')]
        CC <- Y - A - B*X
        C1 <- O[1,1] + O[1,2]*B + O[2,1]*B + O[2,2]*B^2
        C2 <- 2*(O[1,2] + O[2,2]*B)*CC
        C3 <- O[2,2]*CC^2
        K <- -C2/(2*C1)
        SS <- SS + C1*K^2 + C2*K + C3
    }
    SS/2
}

# apply a calibration to a unknown samples
calibrate <- function(lr,fit,dat,PbUstand=NULL,tst=c(337.13,0.18)){
    XY <- lr2XY(lr=lr,dat=dat)
    ns <- length(XY)
    if (is.null(PbUstand))
        Pb6U8stand <- IsoplotR:::age_to_Pb206U238_ratio(tt=tst[1],st=tst[2])
    Ytstand <- log(Pb6U8stand[1])
    sYtstand <- Pb6U8stand[2]/Pb6U8stand[1]
    La <- rep(0,3*ns)
    E <- matrix(0,4*ns+3,4*ns+3)
    E[4*ns+(1:2),4*ns+(1:2)] <- fit$cov
    E[4*ns+3,4*ns+3] <- sYtstand^2
    J <- matrix(0,nrow=3*ns,ncol=4*ns+3)
    for (i in 1:ns){
        i1 <- i      # X
        i2 <- i+ns   # Y
        i3 <- i+2*ns # L7a
        i4 <- i+3*ns # L4a
        o1 <- i      # out: L6a
        o2 <- i+ns   # out: L7a
        o3 <- i+2*ns # out: L4a
        xy <- XY[[i]]
        La[o1] <- xy$x['Y'] + Ytstand - fit$AB['A'] - fit$AB['B']*xy$x['X']
        La[o2] <- xy$x['L7a']
        La[o3] <- xy$x['L4a']
        J[o1,i1] <- -fit$AB['B']   # dL6dX
        J[o1,i2] <- 1              # dL6dY
        J[o1,4*ns+1] <- -1         # dL6dA
        J[o1,4*ns+2] <- -xy$x['X'] # dL6dB
        J[o1,4*ns+3] <- 1          # dL6dYtstand
        J[o2,i3] <- 1              # dL7dL7
        J[o3,i4] <- 1              # dL4dL4
        E[c(i1,i2,i3,i4),c(i1,i2,i3,i4)] <- xy$cov
    }
    out <- list()
    out$snames <- names(lr)
    out$labels <- c('L6a','L7a','L4a')
    out$x <- La
    out$cov <- J %*% E %*% t(J)
    out
}
