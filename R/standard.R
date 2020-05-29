#' @title define an isotopic reference standards
#' @description specify the isotopic composition of a reference
#'     material for SIMS data calibration.
#' @param preset (optional) text string. One of either
#'     \code{'Plesovice'}, \code{'44069'}, \code{'Temora'}, or
#'     \code{'NBS28'}.
#' @param prefix text string marking the first characters of the
#'     sample names that mark the standard analyses
#' @param tst (optional) two-element vector with the age and standard
#'     error of the (presumed concordant) age standard and its
#'     analytical uncertainty.
#' @param val (optional) true isotopic composition of
#'     \eqn{^{206}}Pb/\eqn{^{238}}U- and
#'     \eqn{^{208}}Pb/\eqn{^{232}}Th-ratios of the U-Pb age standard,
#'     or true \eqn{\delta^{18}O} and \eqn{\delta^{16}O} values of the
#'     oxygen isotope reference material.
#' @param cov (optional) the covariance matrix of \code{val}
#' @param common (optional) the common-Pb composition
#'     (\eqn{^{206}}Pb/\eqn{^{204}}Pb and
#'     \eqn{^{208}}Pb/\eqn{^{204}}Pb).
#' @return an object of class \code{standard}
#' @examples
#' data(Cameca,package="simplex")
#' dc <- drift(x=Cameca)
#' lr <- logratios(x=dc)
#' st <- standard(preset="Plesovice")
#' cal <- calibration(lr=lr,stand=st)
#' plot(cal,option=3)
#' @export
standard <- function(preset,prefix=preset,tst,
                     val,cov=matrix(0,length(val),length(val)),common){
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
            if (missing(common)) out$common <- val*0
            else out$common <- common
            out$fetch <- function(dat){
                lrstand(dat)
            }
        }
        out$prefix <- prefix
        class(out) <- 'standard'
    } else if (preset=='Plesovice'){
        out <- standard(tst=c(337.13,0.18),prefix=prefix)
    } else if (preset=='44069'){
        out <- standard(tst=c(424.86,0.25),prefix=prefix)
    } else if (preset=='Temora'){
        out <- standard(tst=c(416.75,0.12),prefix=prefix)
    } else if (preset=='NBS28'){
        out <- standard(val=c(4.79,9.56),cov=diag(c(0.05,0.11))^2,prefix=prefix)
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
        out$common <- dat$stand$common
        names(out$common) <- labels
        out$lr <- log(val)
        J <- diag(1/val)
    } else if (type=="oxygen"){
        if (length(val)==1){ # if d17O is missing
            val <- c(val/2,val)
            J <- matrix(c(0.5,1),2,1)
            cov <- J %*% cov %*% t(J)
        }
        labels <- c("O17O16","O18O16")
        out$ref <- VSMOW()$lr
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
    if (type%in%c('U-Pb','U-Th-Pb')){
        r <- IsoplotR:::age_to_cottle_ratios(tt=tst[1],st=tst[2])
    } else {
        stop('Invalid type argument supplied to age2lr')
    }
    out <- list()
    out$lr <- log(r$x)
    J <- diag(1/r$x)
    rownames(J) <- names(r$x)
    out$cov <- J %*% r$cov %*% t(J)
    out$common <- IsoplotR:::stacey.kramers(tst[1])[,c('i64','i84')]
    names(out$common) <- c('Pb206Pb204','Pb208Pb204')
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
