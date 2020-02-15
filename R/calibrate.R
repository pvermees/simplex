#' @title apply a calibration to sample data
#' @description fit a straight line to log[Pb/U] vs. log[UOx/U] sample
#'     data
#' @details Regresses a line through
#'     X=\eqn{^{238}}U\eqn{^{16}}O\eqn{_x}/\eqn{^{238}}U and
#'     Y=(\eqn{^{206}}Pb-\eqn{^{204}}Pb[\eqn{^{206}}Pb/
#'     \eqn{^{204}}Pb]\eqn{_{\circ}})/\eqn{^{238}}U using the same
#'     slope as previously determined on an age standard
#' @param dat a dataset of class \code{simplex}
#' @param fit the output of \code{calibration}
#' @param syserr propagate the systematic error associated with the
#'     standard (age uncertainty and calibration fit)?
#' @return an \code{IsoplotR} object of class \code{UPb}
#'     (\code{format=5})
#' @examples
#' data(Cameca,package="simplex")
#' stand <- standards(dat=Cameca,prefix='Plesovice',tst=c(337.13,0.18))
#' fit <- calibration(stand,oxide='UO2')
#' unk <- unknowns(dat=Cameca,prefix='Plesovice',invert=TRUE)
#' cal <- calibrate(unk,fit)
#' @export
calibrate <- function(dat,fit,syserr=FALSE){
    # 1. calibrate, calculate (co)variances and partial derivatives
    snames <- names(dat)
    ns <- length(snames)
    X <- matrix(0,ns,8)
    rownames(X) <- snames
    for (sname in snames){
        p <- pars(spot=dat[[sname]],oxide=fit$oxide)
        bg <- get_bg(spot=spot,oxide=fit$oxide)
        X[sname,] <- calibrate_spot(B=fit$AB['b'],p=p,bg=bg)
    }
    # 2. Collate calibrated data into one big covariance structure
    out <- list()
    class(out) <- 'UPb'
    out$format <- 9
    out$names <- rownames(X)
    out$logratios <- c('U238Pb206','Pb207Pb206','Pb204Pb206')
    ns <- length(out$names)
    nr <- length(out$logratios)
    out$x <- rep(0,ns*nr)
    covmat <- matrix(0,ns*nr,ns*nr)
    J <- matrix(0,nrow=ns*nr,ncol=3)
    colnames(J) <- c('A','B','U8Pb6s')
    U8Pb6s <- -log(fit$PbU[1])
    varU8Pb6s <- fit$PbU[2]/fit$PbU[1]
    for (i in 1:ns){
        j <- (i-1)*nr+(1:nr)
        out$x[j[1]] <- U8Pb6s + fit$AB['A'] - X[i,'Ax']
        out$x[j[2:3]] <- X[i,c('Pb76','Pb46')]
        J[j[1],'A'] <- 1
        J[j[1],'B'] <- -X[i,'dAxdBs']
        J[j[1],'U8Pb6s'] <- 1
        covmat[j[1],j[1]] <- X[i,'varAx']
        covmat[j[2],j[2]] <- X[i,'varPb76']
        covmat[j[3],j[3]] <- X[i,'varPb46']
        covmat[j[2],j[3]] <- X[i,'covPb46Pb76']
        covmat[j[3],j[2]] <- covmat[j[2],j[3]]
    }
    E <- matrix(0,3,3)
    E[1:2,1:2] <- fit$cov
    E[3,3] <- varU8Pb6s
    if (syserr) out$cov <- covmat + J %*% E %*% t(J)
    else out$cov <- covmat
    out
}

calibrate_spot <- function(B,p,b0g){
    out <- rep(0,8)
    names(out) <- c('Ax','Pb76','Pb46',
                    'varAx','varPb76','varPb46','covPb46Pb76',
                    'dAxdB')
    init <- log(sum(p$Pb206$c)) - log(sum(p$U238$c))
    out['Ax'] <- stats::optimise(misfit_A,interval=c(-10,10),
                                 B=B,p=p,b0g=b0g)$minimum
    H <- stats::optimHess(par=out['Ax'],fn=misfit_A,B=B,p=p,b0g=b0g)
    out['varAx'] <- solve(H)
    out['dAxdB'] <- misfit_A(A=out['Ax'],B=B,p=p,b0g=b0g,c64=0,deriv=TRUE)
    out['Pb76'] <- getPbLogRatio(p=p,b0g=b0g,num='Pb207',den='Pb206')
    out['Pb46'] <- getPbLogRatio(p=p,b0g=b0g,num='Pb204',den='Pb206')
    HPb <- matrix(0,2,2)
    #HPb[1,1] <- TODO
    #HPb[1,2] <- TODO
    #HPb[2,2] <- TODO
    #HPb[2,1] <- TODO
    covmat <- HPb#solve(-HPb) TODO
    out['varPb46'] <- covmat[1,1]
    out['varPb76'] <- covmat[2,2]
    out['covPb46Pb76'] <- covmat[1,2]
    out
}
