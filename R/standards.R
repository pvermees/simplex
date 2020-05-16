#' @title define the standards in a dataset
#' @description select a subset of standards from a \code{simplex}
#'     dataset. The true isotopic composition of the standard can be
#'     supplied either explicitly, or via its age.
#' @param dat an object of class \code{simplex}
#' @param prefix text string to match
#' @param invert logical.  If \code{TRUE} return samples whose names
#'     do _not_ match
#' @param ratios (optional) true \eqn{^{206}}Pb/\eqn{^{238}}U- and
#'     \eqn{^{208}}Pb/\eqn{^{232}}Th-ratios of the age standard
#' @param ratios.cov (optional) the covariance matrix of \code{DP}
#' @param tst (optional) two-element vector with the age and standard
#'     error of the (presumed concordant) age standard and its
#'     analytical uncertainty
#' @return an object of class \code{standard}
#' @examples
#' data(Cameca,package="simplex")
#' stand <- standards(dat=Cameca,prefix='Plesovice',tst=c(337.13,0.18))
#' @export
standards <- function(dat,prefix,invert=FALSE,type='U-Pb',
                      val,cov=matrix(0,length(val),length(val)),tst){
    out <- list()
    if (missing(val)){
        if (missing(tst)){
            warning('No standard age or composition was supplied.')
            tst <- dat2age(dat,type=type)
        }
        out <- age2lr(tst=tst,type=type)
    } else {
        out <- lrstand(val=val,cov=cov,type=type)
    }
    out$x <- subset_samples(dat=dat,prefix=prefix,invert=invert)
    class(out) <- 'standards'
    out
}
lrstand <- function(val,cov=matrix(0,length(val),length(val)),type="U-Pb"){
    out <- list()
    if (type=="U-Pb"){
        labels <- c("Pb206U238","Pb208Th232")
        out$lr <- log(val)
        J <- diag(1/val)
    } else if (type=="d18O"){
        if (length(val)==1){ # if d17O is missing
            val <- c(val,val/2)
            J <- matrix(c(1,0.5),2,1)
            cov <- J %*% cov %*% t(J)
        }
        labels <- c("O18O16","O17O16")
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
dat2age <- function(dat,type='U-Pb'){
    if (type=='U-Pb'){
        out <- Pb76_to_age(dat)
    } else {
        stop('Invalid type argument supplied to dat2age')
    }
    out
}
age2lr <- function(tst,type='U-Pb'){
    if (type=='U-Pb'){
        r <- IsoplotR:::age_to_cottle_ratios(tt=tst[1],st=tst[2])
    } else {
        stop('Invalid type argument supplied to age2ratios')
    }
    out <- list()
    out$lr <- log(r$x)
    J <- diag(1/r$x)
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



