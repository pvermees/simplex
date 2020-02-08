#' @title plot the raw mass spectrometer data
#' @description Shows the raw data for a single spot in a SIMS
#'     dataset.
#' @param samp one item of a \code{simplex} list object
#' @param fit the output of \code{calibration}
#' @param ... optional parameters to the generic \code{plot} function
#' @return a multi-panel plot
#' @examples
#' data(Cameca,package="simplex")
#' stand <- standard(dat=Cameca,prefix='Plesovice',tst=c(337.13,0.18))
#' samp <- unknown(dat=Cameca,prefix='Qinghu')
#' cal <- calibration(stand,oxide='UO2')
#' plot_timeresolved(samp=samp[[1]],fit=cal)
#' @export
plot_timeresolved <- function(samp,fit=NULL,...){
    if ('unknown' %in% class(samp) & !is.null(fit)){
        calspot <- calibrate_spot(samp,fit=fit)
        cal <- fit
        cal$AB['A'] <- calspot['Ax']
        E <- matrix(0,3,3)
        E[1,1] <- calspot['varAx']
        E[2:3,2:3] <- fit$cov
        J <- matrix(0,2,3)
        J[1,1] <- 1             # dAxdAx
        J[1,3] <- calspot['dAxdBs'] # dAxdBs
        J[2,3] <- 1             # dBsdBs
        cal$cov <- J %*% E %*% t(J)
    } else {
        cal <- fit
    }
    ions <- names(samp$dwelltime)
    np <- length(ions)      # number of plot panels
    nr <- ceiling(sqrt(np)) # number of rows
    nc <- ceiling(np/nr)    # number of columns
    oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    simplex <- NULL
    if (!is.null(fit)){
        simplex <- c('Pb204','Pb206','Pb207','U238',cal$oxide)
        X <- samp$time[,simplex]
        Y <- predict_counts(samp,fit=cal)[,simplex]
    }
    for (ion in ions){
        graphics::plot(samp$time[,ion],samp$counts[,ion],
                       type='p',xlab='',ylab='',...)
        if (!is.null(fit) & ion%in%simplex){
            graphics::lines(X[,ion],Y[,ion])
        }
        graphics::mtext(side=1,text='t',line=2)
        graphics::mtext(side=2,text=ion,line=2)
    }
    graphics::par(oldpar)
}

#' @title calibration plot
#' @description plot Y=(\eqn{^{206}}Pb-\eqn{^{204}}Pb[\eqn{^{206}}Pb/
#'     \eqn{^{204}}Pb]\eqn{_{\circ}})/\eqn{^{238}}U versus
#'     X=\eqn{^{238}}U\eqn{^{16}}O\eqn{_x}/\eqn{^{238}}U
#' @param stand an object of class \code{standard}
#' @param fit the output of \code{calibration}
#' @param labels label the spots with their name (\code{labels=1}) or
#'     number (\code{labels=2}). The default is to use no labels.
#' @return a plot of Y=(\eqn{^{206}}Pb-\eqn{^{204}}Pb[\eqn{^{206}}Pb/
#'     \eqn{^{204}}Pb]\eqn{_{\circ}})/\eqn{^{238}}U versus
#'     X=\eqn{^{238}}U\eqn{^{16}}O\eqn{_x}/\eqn{^{238}}U
#' @examples
#' data(Cameca,package="simplex")
#' Ples <- standards(dat=Cameca,prefix='Plesovice',tst=c(337.13,0.18))
#' cal <- calibration(stand=Ples,oxide='UO2')
#' calplot(stand=Ples,fit=cal)
#' @export
calplot <- function(stand,fit,labels=0,omit=NULL){
    dat <- stand
    if (!is.null(omit)) dat$x <- stand$x[-omit]
    bg <- get_bg(dat,oxide=fit$oxide)
    snames <- names(stand$x)
    nc <- length(snames)
    nr <- nrow(stand$x[[1]]$counts)
    X <- matrix(0,nr,nc)
    Y <- matrix(0,nr,nc)
    colnames(X) <- snames
    colnames(Y) <- snames
    for (sname in snames){
        spot <- stand$x[[sname]]
        p <- pars(spot,oxide=fit$oxide)
        cc <- get_cal_components(p=p,bg=bg,sname=sname)
        b4corr <- log(1 - exp(cc$bdc46)*stand$c64)
        X[,sname] <- cc$bmOU - cc$dcOU + cc$bcO - cc$bcU
        Y[,sname] <- log(p$Pb206$c) - log(p$U238$c) -
            cc$dc6U + cc$bc6 - cc$bcU + b4corr
    }
    tit <- paste0('Y = ',signif(fit$AB['A'],3),'+',
                  signif(fit$AB['B'],3),'X')
    graphics::plot(X,Y,type='n',
                   xlab=paste0('log[',fit$oxide,'/U]'),
                   ylab=paste0('log[Pb','/U]'),main=tit)
    graphics::matlines(X,Y,lty=1,col='grey')
    bg <- rep('black',nc)
    bg[fit$omit] <- 'grey'
    graphics::points(X[1,],Y[1,],pch=21,bg=bg)
    if (labels==1) graphics::text(X[1,],Y[1,],labels=snames,pos=1,col=bg)
    if (labels==2) graphics::text(X[1,],Y[1,],labels=1:length(snames),pos=1,col=bg)
    graphics::points(X[nr,],Y[nr,],pch=21,bg='white')
    xlim <- graphics::par('usr')[1:2]
    graphics::lines(xlim,fit$AB['A']+fit$AB['B']*xlim)
}

predict_counts <- function(samp,fit,c64=0){
    p <- pars(samp,oxide=fit$oxide)
    ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@52@"]]));##:ess-bp-end:##
    bg <- get_bg(samp,oxide=fit$oxide)
    n6 <- p$nU*exp(log_c6U)*p$d6/p$dU
    n7 <- n6*exp(b7+g$Pb*(p$t7-p$t6))*p$d7/p$d6
    n4 <- n6*exp(b4+g$Pb*(p$t4-p$t6))*p$d4/p$d6
    out <- cbind(n4,n6,n7,p$nU,p$nO)
    colnames(out) <- c('Pb204','Pb206','Pb207','U238',fit$oxide)
    out
}
