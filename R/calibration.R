#' @title create a calibration curve
#' @description fit a straight line to log[Pb/U] vs. log[UOx/U] data
#' @details Regresses a line through
#'     X=\eqn{^{238}}U\eqn{^{16}}O\eqn{_x}/\eqn{^{238}}U and
#'     Y=(\eqn{^{206}}Pb-\eqn{^{204}}Pb[\eqn{^{206}}Pb/
#'        \eqn{^{204}}Pb]\eqn{_{\circ}})/\eqn{^{238}}U
#' @param stand a dataset of class \code{simplex}
#' @param oxide either \code{'UO'} or \code{'UO2'}
#' @param c64 the \eqn{^{206}}Pb/\eqn{^{204}}Pb-ratio of the common Pb
#' @param PbU (optional) true \eqn{^{206}}Pb/\eqn{^{238}}U-ratio of
#'     the age standard
#' @param tst (optional) two-element vector with the age and standard
#'     error of the age standard
#' @return a list with the following items:
#' 
#' \code{AB} a vector with the intercept (\code{'A'}) and slope
#' (\code{'B'}) of the power law
#'
#' \code{cov} the covariance matrix of \code{AB}
#'
#' @examples
#' data(Cameca,package="simplex")
#' stand <- standards(dat=Cameca,prefix='Plesovice')
#' fit <- calibration(stand,tst=c(337.13,0.18))
#' @export
calibration <- function(stand,oxide='UO2',c64=18.7,PbU=NULL,tst=NULL){
    fit <- stats::optim(c(A=0,B=1),AB_misfit,dat=stand,oxide=oxide,c64=c64)
    hess <- stats::optimHess(fit$par,AB_misfit,dat=stand,oxide=oxide,c64=c64)
    out <- list()
    out$AB <- fit$par
    out$cov <- solve(hess)
    out$oxide <- oxide
    out$c64 <- c64
    if (is.null(PbU)){
        if (is.null(tst)){
            Pb76 <- get_Pb76(stand)
            warning('No standard age was supplied')
            tst <- IsoplotR:::get.Pb207Pb206.age.default(x=Pb76[1],sx=Pb76[2])
        }
        out$PbU <- IsoplotR:::age_to_Pb206U238_ratio(tt=tst[1],st=tst[2])
    } else {
        out$PbU <- PbU
    }
    out
}

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

AB_misfit <- function(AB,dat,oxide='UO2',c64=18.7){
    out <- 0
    snames <- names(dat)
    for (sname in snames){
        spot <- dat[[sname]]
        out <- out - LL_AB(AB,spot,oxide=oxide,c64=c64)
    }
    out
}

LL_AB <- function(AB,spot,oxide='UO2',c64=18.7){
    p <- pars(spot,oxide=oxide)
    g <- get_gamma(B=AB['B'],p=p)
    a <- get_alpha(AB=AB,p=p,g=g,c64=c64)
    # 1. get count ratios
    n6U <- exp(a$a6 + g$Pb*(p$t6-p$tU))*p$d6/p$dU
    # 2. get proportions
    nt <- length(p$tU)
    den <- n6U + 1
    theta <- cbind(n6U,1)/matrix(rep(den,2),ncol=2)
    colnames(theta) <- c('Pb206','U238')
    counts <- spot$counts[,c('Pb206','U238')]
    # 3. compute multinomial log-likelihood
    sum(counts*log(theta))
}

pars <- function(spot,oxide='UO2'){
    out <- list()
    out$t4 <- spot$time[,'Pb204']
    out$t6 <- spot$time[,'Pb206']
    out$t7 <- spot$time[,'Pb207']
    out$tU <- spot$time[,'U238']
    out$tO <- spot$time[,oxide]
    out$n4 <- spot$counts[,'Pb204']
    out$n6 <- spot$counts[,'Pb206']
    out$n7 <- spot$counts[,'Pb207']
    out$nU <- spot$counts[,'U238']
    out$nO <- spot$counts[,oxide]
    out$c4 <- spot$cps[,'Pb204']
    out$c6 <- spot$cps[,'Pb206']
    out$c7 <- spot$cps[,'Pb207']
    out$cU <- spot$cps[,'U238']
    out$cO <- spot$cps[,oxide]
    out$d4 <- spot$dwelltime['Pb204']
    out$d6 <- spot$dwelltime['Pb206']
    out$d7 <- spot$dwelltime['Pb207']
    out$dU <- spot$dwelltime['U238']
    out$dO <- spot$dwelltime[oxide]
    out
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
