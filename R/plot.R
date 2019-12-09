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
        simplex <- c('Pb204','Pb206','U238',cal$oxide)
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
calplot <- function(stand,fit,labels=0){
    snames <- names(stand$x)
    nc <- length(snames)
    nr <- nrow(stand$x[[1]]$counts)
    X <- matrix(0,nr,nc)
    Y <- matrix(0,nr,nc)
    colnames(X) <- snames
    colnames(Y) <- snames
    for (sname in snames){
        p <- pars(stand$x[[sname]],oxide=fit$oxide)
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
    bg <- rep('black',nc)
    bg[fit$omit] <- 'grey'
    graphics::points(X[1,],Y[1,],pch=21,bg=bg)
    if (labels==1) graphics::text(X[1,],Y[1,],labels=snames,pos=1,col=bg)
    if (labels==2) graphics::text(X[1,],Y[1,],labels=1:length(snames),pos=1,col=bg)
    graphics::points(X[nr,],Y[nr,],pch=21,bg='white')
    xlim <- graphics::par('usr')[1:2]
    graphics::lines(xlim,fit$AB['A']+fit$AB['B']*xlim)
}

predict_counts <- function(samp,fit){
    p <- pars(samp,oxide=fit$oxide)
    g <- get_gamma(B=fit$AB['B'],p=p)
    a <- get_alpha(AB=fit$AB,p=p,g=g,c64=fit$c64)
    n6U <- exp(a$a6 + g$Pb*(p$t6-p$tU))*p$d6/p$dU
    n6 <- p$nU*n6U
    n46 <- (1-exp(a$a4))*(p$d4/p$d6)/fit$c64
    n4 <- p$n6*n46
    out <- cbind(n4,n6,p$nU,p$nO)
    colnames(out) <- c('Pb204','Pb206','U238',fit$oxide)
    out
}

# modified version of filled.contour with ".filled.contour" part replaced with "image"
# function. Note that the color palette is a flipped heat.colors rather than cm.colors
image.with.legend <- function (x = seq(1, nrow(z), length.out = nrow(z)), y = seq(1, 
    ncol(z), length.out=nrow(z)), z, xlim = range(x, finite = TRUE), 
    ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
    levels = pretty(zlim, nlevels), nlevels = 20, color.palette = grDevices::heat.colors, 
    col = rev(color.palette(length(levels) - 1)), plot.title, plot.axes,
    key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
    axes = TRUE, frame.plot = axes, ...) {
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq.int(1, nrow(z), length.out = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
        stop("increasing 'x' and 'y' values expected")
    mar.orig <- (par.orig <- graphics::par(c("mar", "las", "mfrow")))$mar
    on.exit(graphics::par(par.orig))
    w <- (3 + mar.orig[2L]) * graphics::par("csi") * 2.54
    graphics::layout(matrix(c(2, 1), ncol = 2L), widths = c(1, graphics::lcm(w)))
    graphics::par(las = las)
    mar <- mar.orig
    mar[4L] <- mar[2L]
    mar[2L] <- 1
    graphics::par(mar = mar)
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = range(levels),
                xaxs = "i", yaxs = "i")
    graphics::rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
    if (missing(key.axes)) {
        if (axes) 
            graphics::axis(4)
    }
    else key.axes
    graphics::box()
    if (!missing(key.title)) 
        key.title
    mar <- mar.orig
    mar[4L] <- 1
    graphics::par(mar = mar)
    graphics::image(x,y,z,col=col,xlab="",ylab="")
    if (missing(plot.axes)) {
        if (axes) {
            graphics::title(main = "", xlab = "", ylab = "")
            graphics::Axis(x, side = 1)
            graphics::Axis(y, side = 2)
        }
    }
    else plot.axes
    if (frame.plot) 
        graphics::box()
    if (missing(plot.title)) 
        graphics::title(...)
    else plot.title
    invisible()
}
