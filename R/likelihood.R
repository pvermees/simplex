# beta is a matrix
LL <- function(spot,beta){
    nlr <- length(beta$num)
    out <- 0
    ntype <- spot$type[beta$num]
    dtype <- spot$type[beta$den]
    Fc <- (ntype%in%'Fc' & dtype%in%'Fc')
    Em <- (ntype%in%'Em' & dtype%in%'Em')
    Mx <- (ntype%in%'Em' & dtype%in%'Fc') | (ntype%in%'Fc' & dtype%in%'Em')
    if (any(Fc)){
        ll <- LL_Faraday(spot=spot,beta=beta,i=which(Fc))
        out <- out + sum(ll)
    }
    if (any(Em)){
        ll <- LL_SEM(spot=spot,beta=beta,i=which(Em))
        out <- out + sum(ll)
    }
    if (any(Mx)){
        ll <- LL_mixed(spot=spot,beta=beta,i=which(Mx))
        out <- out + sum(ll)
    }
    out
}

LL_Faraday <- function(spot,beta,i){
    num <- beta$num[i]
    den <- beta$den[i]
    X <- log(spot$signal[,num]) - log(spot$signal[,den])
    stable <- (nrow(beta$lr)==1)
    single <- (length(i)==1)
    mu <- beta$lr[,i]
    if (stable){
        if (single) E <- var(X)
        else E <- cov(X)
    } else {
        nt <- nrow(spot$signal)
        if (single) E <- sum((X-mu)^2)/(nt-1)
        else E <- NULL # TODO
    }
    mahalanobis(x=X,center=mu,cov=E)
}

LL_SEM <- function(spot,beta,i){
    0
}

LL_mixed <- function(spot,beta,i){
    0
}

stable <- function(beta){
    if (is.matrix(beta)){
        if (nrow(beta)>1) return(FALSE)
        else return(TRUE)
    } else {
        stop('beta is not a matrix')
    }
}
