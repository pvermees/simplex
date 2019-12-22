# get geometric mean Pb207/Pb206 ratio to estimate
# the standard age if not supplied by the user
get_Pb76 <- function(dat){
    snames <- names(dat)
    ns <- length(snames)
    lPb76 <- rep(0,ns)
    for (i in 1:ns){
        p <- pars(dat[[i]])
        lPb76[i] <- log(sum(p$c7)/sum(p$c6))
    }
    lPb76 <- mean(lPb76)
    slPb76 <- stats::sd(lPb76)/sqrt(ns)
    Pb76 <- exp(lPb76)
    sPb76 <- Pb76*slPb76
    c(Pb76,sPb76)
}

get_gamma <- function(B,p){
    out <- list()
    out$U <- stats::glm(p$nU ~ p$tU,
                        family=stats::poisson(link='log'))$coef[2]
    out$O <- stats::glm(p$nO ~ p$tO,
                        family=stats::poisson(link='log'))$coef[2]
    out$Pb <- B*out$O + (1-B)*out$U
    out
}

get_alpha <- function(AB,p,g,c64){
    out <- list()
    out$a4 <- log(1-(p$c4/p$c6)*c64*exp(g$Pb*(p$t6-p$t4)))
    out$aO <- log(p$cO/p$cU) + g$O*(p$tU-p$tO)
    out$a6 <- AB[1] + AB[2]*out$aO - out$a4
    out
}

logratios2ratios <- function(lr){
    out <- list()
    out$x <- exp(lr$x)
    J <- diag(out$x)
    out$cov <- J %*% lr$cov %*% t(J)
    out
}
