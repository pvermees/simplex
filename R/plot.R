#' @title plot the raw mass spectrometer data
#' @description
#' Shows the raw data for a single spot in a SIMS dataset.
#' @param samp one item of a \code{simplex} list object
#' @param fit the output of \code{calibration}
#' @return a multi-panel plot
#' @examples
#' data(Cameca,package="simplex")
#' plot_timeresolved(Cameca[[1]])
#' @export
plot_timeresolved <- function(samp,fit=NULL){
    ions <- names(samp$dwelltime)
    np <- length(ions)      # number of plot panels
    nr <- ceiling(sqrt(np)) # number of rows
    nc <- ceiling(np/nr)    # number of columns
    graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    simplex <- NULL
    if (!is.null(fit)){
        simplex <- c('Pb204','Pb206','U238',fit$oxide)
        X <- samp$time[,simplex]
        Y <- predict_counts(samp,fit=fit)[,simplex]
    }
    for (ion in ions){
        graphics::plot(samp$time[,ion],samp$counts[,ion],
                       type='p',xlab='',ylab='')
        if (!is.null(fit) & ion%in%simplex){
            graphics::lines(X[,ion],Y[,ion])
        }
        graphics::mtext(side=1,text='t',line=2)
        graphics::mtext(side=2,text=ion,line=2)
    }
}

#' @title calibration plot
#' @description
#' plot Y=(\eqn{^{206}}Pb-\eqn{^{204}}Pb[\eqn{^{206}}Pb/
#'         \eqn{^{204}}Pb]\eqn{_{\circ}})/\eqn{^{238}}U
#'     versus X=\eqn{^{238}}U\eqn{^{16}}O\eqn{_x}/\eqn{^{238}}U
#' @param dat an object of class \code{simplex}
#' @param fit the output of \code{calibration}
#' @return
#' a plot of Y=(\eqn{^{206}}Pb-\eqn{^{204}}Pb[\eqn{^{206}}Pb/
#'              \eqn{^{204}}Pb]\eqn{_{\circ}})/\eqn{^{238}}U
#'     versus X=\eqn{^{238}}U\eqn{^{16}}O\eqn{_x}/\eqn{^{238}}U
#' @examples
#' data(Cameca,package="simplex")
#' Ples <- subset_samples(dat=Cameca,prefix='Plesovice')
#' cal <- calibration(stand=Ples,oxide='UO2')
#' calplot(dat=Ples,fit=cal)
#' @export
calplot <- function(dat,fit){
    snames <- names(dat)
    nc <- length(snames)
    nr <- nrow(dat[[1]]$counts)
    X <- matrix(0,nr,nc)
    Y <- matrix(0,nr,nc)
    colnames(X) <- snames
    colnames(Y) <- snames
    for (sname in snames){
        p <- pars(dat[[sname]],oxide=fit$oxide)
        g <- get_gamma(B=fit$AB['B'],p=p)
        a <- get_alpha(AB=fit$AB,p=p,g=g,c64=fit$c64)
        X[,sname] <- log(p$cO/p$cU) + g$O*(p$tU-p$tO)
        Y[,sname] <- log(p$c6/p$cU) + g$Pb*(p$tU-p$t6) + 
            log(1 - (p$c4/p$c6)*g$Pb*(p$t6-p$t4)*fit$c64)
    }
    tit <- paste0('Y = ',signif(fit$AB['A'],3),'+',
                  signif(fit$AB['B'],3),'X')
    graphics::plot(X,Y,type='n',
                   xlab=paste0('log[',fit$oxide,'/U]'),
                   ylab=paste0('log[Pb','/U]'),main=tit)
    graphics::matlines(X,Y,lty=1,col='grey')
    graphics::points(X[1,],Y[1,],pch=21,bg='black')
    graphics::points(X[nr,],Y[nr,],pch=21,bg='white')
    xlim <- graphics::par('usr')[1:2]
    graphics::lines(xlim,fit$AB['A']+fit$AB['B']*xlim)
}

predict_counts <- function(samp,fit){
    p <- pars(samp,oxide=fit$oxide)
    g <- get_gamma(A=fit$AB['B'],p=p)
    a <- get_alpha(AB=fit$AB,p=p,g=g,c64=fit$c64)
    n6U <- exp(a$a6 + g$Pb*(p$t6-p$tU))*p$d6/p$dU
    n6 <- p$nU*n6U
    n46 <- (1-exp(a$a4))*(p$d4/p$d6)/fit$c64
    n4 <- p$n6*n46
    out <- cbind(n4,n6,p$nU,p$nO)
    colnames(out) <- c('Pb204','Pb206','U238',fit$oxide)
    out
}
