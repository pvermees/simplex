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
#'     or the true \eqn{\delta}-values of the stable isotope reference
#'     material.
#' @param cov (optional) the covariance matrix of \code{val}
#' @param common (optional) the common-Pb composition
#'     (\eqn{^{206}}Pb/\eqn{^{204}}Pb and
#'     \eqn{^{208}}Pb/\eqn{^{204}}Pb).
#' @return an object of class \code{standard}
#' @examples
#' data(Cameca_UPb,package="simplex")
#' dc <- drift(x=Cameca_UPb)
#' lr <- logratios(x=dc)
#' st <- standard(preset="Plesovice")
#' cal <- calibration(lr=lr,stand=st)
#' plot(cal,option=3)
#' @export
standard <- function(preset,prefix=preset,tst,val,
                     cov=matrix(0,length(val),length(val)),common){
    if (missing(preset)){
        preset <- 'other'
        out <- list()
        if (missing(val)){
            if (missing(tst)){
                stop("Both the age and the isotopic composition ",
                     "of the sample are missing.")
            } else {
                out$tst <- tst
                out$fetchfun <- "age2lr"
            }
        } else {
            out$val <- val
            out$cov <- cov
            if (!missing(common)) out$common <- common
            out$fetchfun <- "lrstand"
        }
        out$prefix <- prefix
        class(out) <- 'standard'
    } else if (preset=='Plesovice'){
        out <- standard(tst=c('t'=337.13,'s[t]'=0.18),prefix=prefix)
    } else if (preset=='Qinghu'){
        out <- standard(tst=c('t'=159.5,'s[t]'=0.1),prefix=prefix)
    } else if (preset=='44069'){
        out <- standard(tst=c('t'=424.86,'s[t]'=0.25),prefix=prefix)
    } else if (preset=='Temora'){
        out <- standard(tst=c('t'=416.75,'s[t]'=0.12),prefix=prefix)
    } else if (preset=='NBS28'){
        out <- standard(val=c('d17O'=4.79,'d18O'=9.56),
                        cov=diag(c(0.05,0.11))^2,
                        prefix=prefix)
    } else if (preset=='Sonora'){ # temporary value
        out <- standard(val=c('d33S'=0.83,'d34S'=1.61,'d36S'=3.25),
                        cov=diag(c(0.03,0.08,0.03))^2,
                        prefix=prefix)
    } else {
        stop("Invalid input to standard(...).")
    }
    out$name <- preset
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
        out$lr <- log(1 + val/1000) + out$ref
        J <- diag(1/(1000 + val))
    } else if (type=="sulphur"){
        if (length(val)==2){ # if d36S is missing
            val <- c(val,val[2]*2)
            J <- rbind(c(1,0),c(0,1),c(0,2))
            cov <- J %*% cov %*% t(J)
        }
        labels <- c("S33S32","S34S32","S36S32")
        out$ref <- troilite()$lr
        out$lr <- log(1 + val/1000) + out$ref
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
