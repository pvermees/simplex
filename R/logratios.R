logratios <- function(dat,fit){
    Pb6U8 <- calibrate_Pb206U238(dat,fit)
    ns <- length(Pb6U8$snames)
    Pb6U8
}
calibrate_Pb206U238 <- function(dat,fit){
    out <- list()
    out$snames <- names(dat)
    ns <- length(out$snames)
    nt <- nrow(dat[[1]]$time)
    Pb6U8 <- rep(0,ns*nt)
    J <- matrix(-1,ns*nt,2)
    colnames(J) <- c('A','B')
    for (i in 1:ns){
        spot <- dat[[i]]
        LJg <- get_Pb206U238(spot,fit)
        j <- (i-1)*nt + 1:nt
        Pb6U8[j] <- LJg$Pb206U238
        J[j,'B'] <- LJg$dPb206U238dB
    }
    out$x <- Pb6U8
    out$cov <- J %*% fit$cov %*% t(J)
    out
}

get_Pb206U238 <- function(spot,fit){
    p <- pars(spot,oxide=fit$oxide)
    g <- get_gamma(B=fit$AB['B'],p=p)
    out <- list()
    Y <- log(p$c6/p$cU) - g$Pb*(p$t6-p$tU)
    X <- log(p$cO/p$cU) - g$O*(p$tO-p$tU)
    out <- list()
    out$Pb206U238 <- log(fit$PbU[1]) + Y - fit$AB['A'] - fit$AB['B']*X
    out$dPb206U238dB <- -X
    out
}

