#' @title define the standards in a dataset
#' @description select a subset of standards from a \code{simplex}
#'     dataset. The true isotopic composition of the standard can be
#'     supplied either explicitly, or via its age.
#' @param dat an object of class \code{simplex}
#' @param val (optional) true \eqn{^{206}}Pb/\eqn{^{238}}U- and
#'     \eqn{^{208}}Pb/\eqn{^{232}}Th-ratios of the age standard
#' @param cov (optional) the covariance matrix of \code{DP}
#' @param tst (optional) two-element vector with the age and standard
#'     error of the (presumed concordant) age standard and its
#'     analytical uncertainty
#' @return an object of class \code{standard}
#' @examples
#' data(Cameca,package="simplex")
#' stand <- standards(dat=Cameca,prefix='Plesovice',tst=c(337.13,0.18))
#' @export
standard <- function(preset,tst,val,cov=matrix(0,length(val),length(val))){
    if (missing(preset)){
        out <- list()
        if (missing(val)){
            if (missing(tst)){
                out$fetch <- function(dat){
                    temp <- dat
                    temp$tst <- dat2age(dat)
                    age2lr(temp)
                }
            } else {
                out$tst <- tst
                out$fetch <- function(dat){
                    age2lr(dat)
                }
            }
        } else {
            out$val <- val
            out$cov <- cov
            out$fetch <- function(dat){
                lrstand(dat)
            }
        }
        class(out) <- 'standard'
    } else if (preset=='Plesovice'){
        out <- standard(tst=c(337.13,0.18))
    } else if (preset=='44069'){
        out <- standard(tst=c(424.86,0.25))
    } else if (preset=='Temora'){
        out <- standard(tst=c(416.75,0.12))
    } else if (preset=='NBS28'){
        out <- standard(val=c(4.79,9.56),cov=diag(c(0.05,0.11))^2)
    } else {
        stop("Invalid input to standard(...).")
    }
    out
}
lrstand <- function(dat){
    val <- dat$stand$val
    cov <- dat$stand$cov
    type <- datatype(dat)
    out <- list()
    if (type=="U-Pb"){
        labels <- c("Pb206U238","Pb208Th232")
        out$lr <- log(val)
        J <- diag(1/val)
    } else if (type=="oxygen"){
        if (length(val)==1){ # if d17O is missing
            val <- c(val/2,val)
            J <- matrix(c(0.5,1),2,1)
            cov <- J %*% cov %*% t(J)
        }
        labels <- c("O17O16","O18O16")
        out$lr <- log(1 + val/1000) + VSMOW()$lr
        J <- diag(1/(1000 + val))
    } else {
        stop('Invalid type argument supplied to lrstand')
    }
    out$cov <- J %*% cov %*% t(J)
    names(out$lr) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    out
}
dat2age <- function(dat){
    if (type=='U-Pb'){
        out <- Pb76_to_age(dat)
    } else {
        stop('Invalid type argument supplied to dat2age')
    }
    out
}
age2lr <- function(dat){
    type <- datatype(dat)
    tst <- dat$stand$tst
    if (type=='U-Pb'){
        r <- IsoplotR:::age_to_cottle_ratios(tt=tst[1],st=tst[2])
    } else {
        stop('Invalid type argument supplied to age2lr')
    }
    out <- list()
    out$lr <- log(r$x)
    J <- diag(1/r$x)
    rownames(J) <- names(r$x)
    out$cov <- J %*% r$cov %*% t(J)
    out
}
# get geometric mean Pb207/Pb206 ratio to estimate
# the standard age if not supplied by the user
Pb76_to_age <- function(dat){
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
    IsoplotR:::get.Pb207Pb206.age.default(x=Pb76,sx=sPb76)
}



