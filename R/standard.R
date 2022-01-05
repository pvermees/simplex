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
standard <- function(preset,tst,measured,del,ref){
    if (!missing(preset)){
        if (preset=='Plesovice'){
            geochron <- TRUE
            out <- age2stand(tst=c(337.13,0.18))
        } else if (preset=='Qinghu'){
            geochron <- TRUE
            out <- age2stand(tst=c(159.5,0.1))
        } else if (preset=='44069'){
            geochron <- TRUE
            out <- age2stand(tst=c(424.86,0.25))
        } else if (preset=='Temora'){
            geochron <- TRUE
            out <- age2stand(tst=c(416.75,0.12))
        } else if (preset=='NBS28'){
            geochron <- FALSE
            del <- delhelper(num=c('O17','O18'),
                             den='O16',
                             val=c(4.79,9.56),
                             cov=diag(c(0.05,0.11))^2)
            out <- del2stand(del,ref=VSMOW())
        } else if (preset=='Sonora'){ # temporary value
            geochron <- FALSE
            del <- delhelper(num=c('S33','S34','S36'),
                             den='S32',
                             val=c(0.83,1.61,3.25),
                             cov=diag(c(0.03,0.08,0.03))^2)
            out <- del2stand(del,ref=troilite())
        } else {
            stop("Invalid input to standard(...).")
        }
        out$preset <- preset
    } else if (!missing(tst)){
        geochron <- TRUE
        out <- age2stand(tst)
    } else if (!missing(del) & !missing(ref)){
        geochron <- FALSE
        out <- del2stand(del,ref)
    } else {
        stop('Illegal input to standard().')
    }
    if (geochron){
        if (missing(measured)) measured <- FALSE
        out$measured <- measured
        if (out$measured) {
            out$val <- out$val[-c(1:2)]
            out$cov <- out$cov[-c(1:2),-c(1:2)]
            if (!missing(preset)){
                if (preset=='Plesovice'){
                    out$val['Pb206/U238'] <- log(0.053707)
                    out$cov['Pb206/U238','Pb206/U238'] <- 0.02
                    out$cov['Pb206/U238','Pb208/Th232'] <- 0
                    out$cov['Pb208/Th232','Pb206/U238'] <- 0
                }
            }
        }
    } else {
        out$measured <- TRUE
    }
    out
}

delhelper <- function(num,den,val,cov){
    ratios <- paste0(num,'/',den)
    out <- list(val=val,cov=cov)
    names(out$val) <- ratios
    rownames(out$cov) <- ratios
    colnames(out$cov) <- ratios
    out
}

age2stand <- function(tst){
    num <- c('Pb204','Pb204','Pb206','Pb208')
    den <- c('Pb206','Pb208','U238','Th232')
    ratios <- paste0(num,'/',den)
    common <- IsoplotR:::stacey.kramers(tst[1])
    D <- IsoplotR::mclean(tt=tst[1])
    val <- c(-log(common[,c('i64','i84')]),log(D$Pb206U238),log(D$Pb208Th232))
    J <- rbind(0,0,D$dPb206U238dt/D$Pb206U238,D$dPb208Th232dt/D$Pb208Th232)
    covmat <- J%*%(tst[2]^2)%*%t(J)
    names(val) <- ratios
    rownames(covmat) <- ratios
    colnames(covmat) <- ratios
    list(tst=tst,val=val,cov=covmat,measured=FALSE)
}

del2stand <- function(del,ref){
    keep <- (names(ref$val) %in% names(del$val))
    if (!any(keep)) stop('Standard isotopes must match reference.')
    out <- list()
    out$del <- list(val=del$val,cov=del$cov)
    out$ref <- list(preset=ref$preset,val=ref$val,cov=ref$cov)
    out$val <- log(1 + del$val/1000) + ref$val[keep]
    J <- diag(sum(keep))/(1000 + del$val)
    rownames(J) <- names(del$val)
    out$cov <- J %*% data.matrix(del$cov) %*% t(J)
    out
}

VSMOW <- function(){
    O678 <- c(0.3799e-3,2.00520e-3)
    relerr <- c(1.6e-3,0.43e-3)/c(0.3799,2.00520)
    out <- delhelper(num = c('O17','O18'),den = 'O16',
                     val=log(O678),cov=diag(relerr^2))
    out$preset <- 'VSMOW'
    out
}

troilite <- function(){
    S2346 <- c(126.948,22.6436,6515)
    relerr <- c(0.047,0.0020,20)/S2346
    out <- delhelper(num=c('S33','S34','S36'),den='S32',
                     val=-log(S2346),cov=diag(relerr^2))
    out$preset <- 'troilite'
    out
}

skeletonstand <- function(lr,measured=TRUE){
    num <- lr$method$num
    den <- lr$method$den
    geochron <- ( all(c("Pb206","U238") %in% c(num,den)) |
                  all(c("Pb208","Th232") %in% c(num,den)) )
    num_is_element <- (element(num) %in% elements())
    den_is_element <- (element(den) %in% elements())
    ratios <- paste0(num,'/',den)
    if (geochron & measured){
        hasPb204 <- (num %in% "Pb204") | (den %in% "Pb204")
        selection <- ratios[num_is_element & den_is_element & !hasPb204]
    } else {
        selection <- ratios[num_is_element & den_is_element]
    }
    n2r <- length(selection)
    out <- list()
    out$measured <- measured
    out$val <- rep(0,n2r)
    out$cov <- matrix(0,n2r,n2r)
    names(out$val) <- selection
    rownames(out$cov) <- selection
    colnames(out$cov) <- selection
    out
}
