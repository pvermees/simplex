#' @title create a calibration curve
#' @description fit a straight line to log[Pb/U] vs. log[UOx/U] data
#' @details Regresses a line through
#'     X=\eqn{^{238}}U\eqn{^{16}}O\eqn{_x}/\eqn{^{238}}U and
#'     Y=(\eqn{^{206}}Pb-\eqn{^{204}}Pb[\eqn{^{206}}Pb/
#'     \eqn{^{204}}Pb]\eqn{_{\circ}})/\eqn{^{238}}U
#' @param stand a dataset of class \code{simplex}
#' @param oxide either \code{'UO'} or \code{'UO2'}
#' @param omit indices of measurements to be omitted from the
#'     calibration
#' @return a list with the following items:
#' 
#' \code{AB} a vector with the intercept (\code{'A'}) and slope
#' (\code{'B'}) of the power law
#'
#' \code{cov} the covariance matrix of \code{AB}
#'
#' @examples
#' data(Cameca,package="simplex")
#' stand <- standards(dat=Cameca,prefix='Plesovice',tst=c(337.13,0.18))
#' fit <- calibration(stand,oxide='UO2',omit=2)
#' calplot(stand,fit)
#' @export
calibration <- function(stand,oxide='UO2',omit=NULL){
    dat <- stand
    if (!is.null(omit)) dat$x <- stand$x[-omit]
    B0G <- get_B0G(dat,oxide='UO2')
    init <- initAB(stand=stand,oxide=oxide)
    fit <- stats::optim(init,misfit_AB,stand=stand,oxide=oxide,B0G=B0G)
    hess <- stats::optimHess(fit$par,misfit_AB,stand=stand,oxide=oxide,B0G=B0G)
    out <- list()
    out$omit <- omit
    out$AB <- fit$par
    out$cov <- solve(hess)
    out$oxide <- oxide
    out$c64 <- stand$c64
    out$PbU <- stand$PbU
    out
}

initAB <- function(stand,oxide='UO2'){
    X <- NULL
    Y <- NULL
    snames <- names(stand$x)
    for (sname in snames){
        p <- pars(spot=stand$x[[sname]],oxide=oxide)
        bOU <- log(p$O$c) - log(p$U238$c)
        b6U <- log(p$Pb206$c) - log(p$U238$c)
        X <- c(X,bOU)
        Y <- c(Y,b6U)
    }
    out <- lm(Y ~ X)$coef
    names(out) <- c('A','B')
    out
}

misfit_AB <- function(AB,stand,oxide='UO2',B0G){
    snames <- names(stand$x)
    out <- 0
    for (sname in snames){
        p <- pars(spot=stand$x[[sname]],oxide=oxide)
        b0g <- B0G[[sname]]
        out <- out + misfit_A(A=AB[1],B=AB[2],p=p,b0g=b0g,c64=stand$c64)
    }
    out
}
misfit_A <- function(A,B,p,b0g,c64=0,deriv=FALSE){
    # predict Pb206/U238 cps logratio (b6Upc):
    bOUmc <- log(p$O$c) - log(p$U238$c)
    bdcorrUO <- A2Corr(p=p,b0g=b0g,num='O',den='U238')
    X <- bOUmc - bdcorrUO
    RHS <- A + B * X
    bdcorr46 <- A2Corr(p=p,b0g=b0g,num='Pb204',den='Pb206')
    b46 <- b0g['Pb204','b0'] - b0g['Pb206','b0'] - bdcorr46
    b4corr <- log(1 - exp(b46)*c64)
    bdcorr6U <- A2Corr(p=p,b0g=b0g,num='Pb206',den='U238')
    b6Upc <- RHS + bdcorr6U - b4corr
    # calculate the log-likelihood
    b6Upn <- b6Upc + log(p$Pb206$d) - log(p$U238$d)
    n6 <- p$Pb206$n
    nU <- p$U238$n
    LL <- n6*b6Upn - (n6+nU)*log(1+exp(b6Upn))
    if (deriv){
        db6UpcdA <- 1
        db6UpcdB <- X
        num <- exp(b6Upc)
        den <- 1 + exp(b6Upc)
        dnumdA <- exp(b6Upc)*db6UpcdA
        dnumdB <- exp(b6Upc)*db6UpcdB
        dLdA <- db6UpcdA*(n6-(n6+nU)*num/den)
        d2LdA2 <- -db6UpcdA*dnumdA*(n6+nU)/(den^2)
        d2LdAdB <- -db6UpcdA*dnumdB*(n6+nU)/(den^2)
        out <- -sum(d2LdAdB)/sum(d2LdA2)
    } else {
        out <- -sum(LL)
    }
    out
}
