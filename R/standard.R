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
#' @param val (optional) true isotopic composition (as logratios) of
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
standard <- function(x){
    if (x=='Plesovice'){
        out <- age2stand(tst=c(337.13,0.18))
    } else if (x=='Qinghu'){
        out <- age2stand(tst=c(159.5,0.1))
    } else if (x=='44069'){
        out <- age2stand(tst=c(424.86,0.25))
    } else if (x=='Temora'){
        out <- age2stand(tst=c(416.75,0.12))
    } else if (x=='NBS28'){
        del <- data.frame(num=c('O17','O18'),den='O16',
                          val=c(4.79,9.56),cov=diag(c(0.05,0.11))^2)
        out <- del2stand(del,ref=VSMOW())
    } else if (x=='Sonora'){ # temporary value
        del <- data.frame(num=c('S33','S34','S36'),den='S32',
                          val=c(4.79,9.56),cov=diag(c(0.05,0.11))^2)
        out <- del2stand(ppdel,ref=troilite())
    } else {
        stop("Invalid input to standard(...).")
    }
    out
}

age2stand <- function(tst,pairing=NULL){
    num <- c('Pb204','Pb204','Pb207','Pb206','Pb208')
    den <- c('Pb206','Pb208','Pb206','U238','Th232')
    ratios <- paste0(num,'/',den)
    common <- IsoplotR:::stacey.kramers(tst[1])
    D <- IsoplotR::mclean(tt=tst[1])
    val <- c(-log(common[,c('i64','i84')]),
             D$Pb207Pb206,D$Pb206U238,D$Pb208Th232)
    J <- rbind(0,0,D$dPb207Pb206dt,D$dPb206U238dt,D$dPb208Th232dt)
    data.frame(ratios=ratios,val=val,cov=J%*%(tst[2]^2)%*%t(J))
}

del2stand <- function(del,ref){
    keep <- (ref$num %in% del$num) & (ref$den %in% del$den)
    if (!any(keep)) stop('Standard isotopes must match reference.')
    out <- data.frame(ratios=paste0(del$num,'/',del$den))
    out$val <- log(1 + del$val/1000) + ref$val[keep]
    J <- diag(sum(keep))/(1000 + del$val)
    covmat <- data.matrix(del[,4:ncol(del)])
    out$cov <- J %*% covmat %*% t(J)
    out
}

VSMOW <- function(){
    O678 <- c(0.3799e-3,2.00520e-3)
    relerr <- c(1.6e-3,0.43e-3)/c(0.3799,2.00520)
    data.frame(
        num = c('O17','O18'),den = 'O16',
        val=log(O678),cov=diag(relerr^2))
}

troilite <- function(){
    S2346 <- c(126.948,22.6436,6515)
    relerr <- c(0.047,0.0020,20)/S2346
    data.frame(num=c('S33','S34','S36'),den='S32',
               val=-log(S2346),cov=diag(relerr^2))
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
        out$ref <- 
        out$lr <- log(1 + val/1000) + out$ref
        J <- diag(1/(1000 + val))
    } else if (type=="sulphur"){
        if (length(val)==2){ # if d36S is missing
            val <- c(val,val[2]*2)
            J <- rbind(c(1,0),c(0,1),c(0,2))
            cov <- J %*% cov %*% t(J)
        }
        labels <- c("S33S32","S34S32","S36S32")
        out$ref <- 
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
