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
calplot <- function(dat,fit,cal,labels=0,omit=NULL){
    # 1. set up the plot variables
    if (is(dat,'standard')) d <- dat$x
    else d <- dat
    if (!is.null(omit)) d <- d[-omit]
    snames <- names(d)
    nc <- length(snames)
    nr <- nrow(d[[1]]$counts)
    X <- matrix(0,nr,nc)
    Y <- matrix(0,nr,nc)
    colnames(X) <- snames
    colnames(Y) <- snames
    # 2. calculate plot coordinates:
    for (sname in snames){
        spot <- d[[sname]]
        p <- pars(spot,parent=fit$parent,daughter=fit$daughter,oxide=fit$oxide)
        b0g <- get_b0g(spot=d[[sname]],parent=fit$parent,oxide=fit$oxide)
        XY <- getCalXY(p=p,b0g=b0g,cD4=fit$cD4)
        X[,sname] <- XY$X
        Y[,sname] <- XY$Y
    }
    # 3. create the actual plot
    tit <- paste0('Y = ',signif(fit$AB['A'],3),'+',
                  signif(fit$AB['B'],3),'X')
    graphics::plot(X,Y,type='n',
                   xlab=paste0('log[',fit$oxide,'/',fit$parent,']'),
                   ylab=paste0('log[',fit$daughter,'/',fit$parent,']'),main=tit)
    graphics::matlines(X,Y,lty=1,col='grey')
    bg <- rep('black',nc)
    bg[fit$omit] <- 'grey'
    graphics::points(X[1,],Y[1,],pch=21,bg=bg)
    if (labels==1) graphics::text(X[1,],Y[1,],labels=snames,pos=1,col=bg)
    if (labels==2) graphics::text(X[1,],Y[1,],labels=1:length(snames),pos=1,col=bg)
    graphics::points(X[nr,],Y[nr,],pch=21,bg='white')
    xlim <- graphics::par('usr')[1:2]
    if (missing(cal)){
        graphics::lines(xlim,fit$AB['A']+fit$AB['B']*xlim)
    } else {
        ns <- length(snames)
        for (i in 1:ns){
            sname <- snames[i]
            xlim <- range(X[,sname])
            dp <- paste0(fit$daughter,fit$parent)
            dA <- cal$x[i] - log(fit$DP[dp])
            ylim <- fit$AB['A']+dA+fit$AB['B']*xlim
            graphics::lines(xlim,ylim)
        }
    }
}

#' @title plot the raw mass spectrometer data
#' @description Shows the raw data for a single spot in a SIMS
#'     dataset.
#' @param spot one item of a \code{simplex} list object
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
plot_timeresolved <- function(spot,fit=NULL,...){
    ions <- names(spot$dwelltime)
    np <- length(ions)      # number of plot panels
    nr <- ceiling(sqrt(np)) # number of rows
    nc <- ceiling(np/nr)    # number of columns
    oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    simplex <- NULL
    if (!is.null(fit)){
        predictable <- c('bkg','Pb204','Pb206','Pb207','U238',fit$oxide)
        plottable <- predictable %in% colnames(spot$time)
        simplex <- predictable[plottable]
        X <- spot$time[,simplex]
        Y <- predict_counts(spot,fit=fit)[,simplex]
    }
    for (ion in ions){
        graphics::plot(spot$time[,ion],spot$counts[,ion],
                       type='p',xlab='',ylab='',...)
        if (!is.null(fit) & ion%in%simplex){
            graphics::lines(X[,ion],Y[,ion])
        }
        graphics::mtext(side=1,text='t',line=2)
        graphics::mtext(side=2,text=ion,line=2)
    }
    graphics::par(oldpar)
}

predict_counts <- function(spot,fit){
    p <- pars(spot=spot,parent=fit$parent,oxide=fit$oxide)
    b0g <- get_b0g(spot=spot,oxide=fit$oxide)
    calspot <- calibrate_spot(B=fit$AB['B'],p=p,b0g=b0g)
    b76pc <- calspot['Pb76'] + A2Corr(p=p,b0g=b0g,num='Pb207',den='Pb206')
    b76pn <- b76pc + log(p$Pb207$d) - log(p$Pb206$d)
    n6 <- (p$Pb206$n+p$Pb207$n)/(1+exp(b76pn))
    n7 <- (p$Pb206$n+p$Pb207$n)*exp(b76pn)/(1+exp(b76pn))
    b86pc <- calspot['Pb76'] + A2Corr(p=p,b0g=b0g,num='Pb207',den='Pb206')
    b86pn <- b86pc + log(p$Pb207$d) - log(p$Pb206$d)
    n8 <- (p$Pb206$n+p$Pb207$n)*exp(b76pn)/(1+exp(b76pn))
    b46pc <- calspot['Pb46'] + A2Corr(p=p,b0g=b0g,num='Pb204',den='Pb206')
    b46pn <- b46pc + log(p$Pb204$d) - log(p$Pb206$d)
    n4 <- (p$Pb206$n+p$Pb204$n)*exp(b46pn)/(1+exp(b46pn))
    bkg <- exp(b0g['blank','b0'])*p$blank$d
    out <- cbind(bkg,n4,n6,n7,n8,p[[p$parent]]$n,p[[p$oxide]]$n)
    colnames(out) <- c('bkg','Pb204','Pb206','Pb207','Pb208',p$parent,p$oxide)
    out
}
