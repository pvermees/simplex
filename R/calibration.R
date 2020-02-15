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
    fit <- stats::optim(init,AB_misfit,stand=stand,oxide=oxide,B0G=B0G)
    hess <- stats::optimHess(fit$par,AB_misfit,stand=stand,oxide=oxide,B0G=B0G)
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

AB_misfit <- function(AB,stand,oxide='UO2',B0G){
    snames <- names(stand$x)
    out <- 0
    for (sname in snames){
        p <- pars(spot=stand$x[[sname]],oxide=oxide)
        b0g <- B0G[[sname]]
        # predict Pb206/U238 cps logratio (b6Upc):
        bOUmc <- log(p$O$c) - log(p$U238$c)
        bdcorrUO <- A2Corr(p=p,b0g=b0g,num='O',den='U238')
        X <- bOUmc - bdcorrUO
        RHS <- AB[1] + AB[2] * X
        bdcorr46 <- A2Corr(p=p,b0g=b0g,num='Pb204',den='Pb206')
        b46 <- b0g['Pb204','b0'] - b0g['Pb206','b0'] - bdcorr46
        b4corr <- log(1 - exp(b46)*stand$c64)
        bdcorr6U <- A2Corr(p=p,b0g=b0g,num='Pb206',den='U238')
        b6Upc <- RHS + bdcorr6U - b4corr
        # calculate the log-likelihood
        b6Upn <- b6Upc + log(p$Pb206$d) - log(p$U238$d)
        LL <- LLbinom(bn=b6Upn,nnum=p$Pb206$n,nden=p$U238$n)
        out <- out - sum(LL)
    }
    out
}

# bg = matrix of beta-gamma values for a single spot
get_cal_components <- function(p,b0g){
    out <- list()
    tU <- p$U238$t
    # bm = measured beta (logratio of cps):
    out$bmOU <- log(p$O$c) - log(p$U238$c)
    out$bm76 <- log(p$Pb207$c) - log(p$Pb206$c)
    out$bm46 <- b0g['Pb204','b'] - b0g['Pb206','b']
    # bc = blank correction:
    out$bcO <- blank_correct(b0g=b0g,tt=tU,mass='O')
    out$bcU <- blank_correct(b0g=b0g,tt=tU,mass='U238')
    out$bc4 <- blank_correct(b0g=b0g,tt=tU,mass='Pb204')
    out$bc6 <- blank_correct(b0g=b0g,tt=tU,mass='Pb206')
    out$bc7 <- blank_correct(b0g=b0g,tt=tU,mass='Pb207')
    # dc = drift correction
    out$dcOU <- b0g['O','g']*(p$O$t-p$U238$t)
    out$dc46 <- b0g['Pb206','g']*(p$Pb204$t-p$Pb206$t)
    out$dc6U <- b0g['Pb206','g']*(p$Pb206$t-p$U238$t)
    out$dc76 <- b0g['Pb206','g']*(p$Pb207$t-p$Pb206$t)
    # blank and drift corrected Pb/Pb logratios
    out$bdc46 <- out$bm46 + out$dc46 + out$bc4 - out$bc6
    out$bdc76 <- out$bm76 - out$dc76 + out$bc7 - out$bc6
    out
}
