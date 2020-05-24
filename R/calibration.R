calibration <- function(lr,stand,prefix=NULL,snames=NULL,i=NULL,invert=FALSE,t=0){
    out <- lr
    dat <- subset(x=lr,prefix=prefix,snames=snames,i=i,invert=invert)
    if (stable(lr)) out$cal <- stable_calibration(lr=dat)
    else out$cal <- geochron_calibration(lr=dat,t=t,DC0=stand$DC0)
    out$stand <- stand
    class(out) <- append('calibration',class(out))
    out
}

stable_calibration <- function(lr){
    LL <- function(par,lr){
        out <- 0
        np <- length(par)
        snames <- names(lr$x)
        for (sname in snames){
            X <- lr$x[[sname]]$lr$b0g[1:np]
            E <- lr$x[[sname]]$lr$cov[1:np,1:np]
            LL <- mahalanobis(x=X,center=par,cov=E)
            out <- out + LL
        }
        out
    }
    ni <- length(lr$m$num)
    init <- lr$x[[1]]$lr$b0g[1:ni]
    wtdmean <- optim(init,fn=LL,gr=NULL,method='BFGS',hessian=TRUE,lr=lr)
    out <- list()
    out$snames <- names(lr$x)
    out$x <- wtdmean$par
    out$cov <- solve(wtdmean$hessian)
    out
}

geochron_calibration <- function(lr,t=0,DC0,...){
    out <- list()
    oxide <- lr$m$oxide
    type <- datatype(lr)
    if (type=="U-Pb"){
        out[[type]] <- geocal(lr,oxide=oxide['U'],t=t,
                              type=type,DC0=DC0['Pb206Pb204'])
    } else if (datatype(lr)=="U-Th-Pb"){
        out[['U-Pb']] <- geocal(lr,oxide=oxide['U'],t=t,
                                type="U-Pb",DC0=DC0['Pb206Pb204'])
        out[['Th-Pb']] <- geocal(lr,oxide=oxide['Th'],t=t,
                                 type="Th-Pb",DC0=DC0['Pb208Pb204'])
    } else {
        stop("Invalid data type.")
    }
    out
}

geocal <- function(lr,oxide,t,type,DC0){
    if (type=='U-Pb'){
        num <- c(oxide,'Pb206','Pb204')
        den <- c('U238','U238','Pb206')
    } else if (type=='Th-Pb'){
        num <- c(oxide,'Pb208','Pb204')
        den <- c('Th232','Th232','Pb208')
    } else {
        stop("Invalid data type.")
    }
    out <- list(num=num,den=den,oxide=oxide,t=t)
    out$DC0 <- DC0
    out$york <- beta2york(lr=lr,t=t,num=num,den=den,DC0=DC0)
    out$fit <- IsoplotR:::york(out$york)
    out
}

beta2york <- function(lr,t=0,num=c('UO2','Pb206','Pb204'),
                      den=c('U238','U238','Pb206'),
                      snames=NULL,i=NULL,DC0=0){
    if (is.null(snames)) snames <- names(lr$x)
    if (!is.null(i)) snames <- snames[i]
    ns <- length(snames)
    out <- matrix(0,nrow=ns,ncol=5)
    colnames(out) <- c('X','sX','Y','sY','rXY')
    rownames(out) <- snames
    for (sname in snames){
        LR <- lr$x[[sname]]$lr
        b0g <- LR$b0g
        b0gnames <- names(b0g)
        nb0g <- length(b0gnames)
        OP <- paste0(num[1],'/',den[1])
        DP <- paste0(num[2],'/',den[2])
        ib0OP <- which(b0gnames %in% paste0('b0[',OP,']'))
        ib0DP <- which(b0gnames %in% paste0('b0[',DP,']'))
        igOP <- which(b0gnames %in% paste0('g[',OP,']'))
        igDP <- which(b0gnames %in% paste0('g[',DP,']'))
        bOP <- b0g[ib0OP] + b0g[igOP]*hours(t)
        bDP <- b0g[ib0DP] + b0g[igDP]*hours(t)
        X <- bOP
        Y <- bDP
        J <- matrix(0,2,length(LR$b0g))
        J[1,ib0OP] <- 1
        J[1,igOP] <- hours(t)
        J[2,ib0DP] <- 1
        J[2,igDP] <- hours(t)
        if (length(num)>2){
            CD <- paste0(num[3],'/',den[3])
            ib0CD <- which(b0gnames %in% paste0('b0[',CD,']'))
            igCD <- which(b0gnames %in% paste0('g[',CD,']'))
            bCD <- b0g[ib0CD] + b0g[igCD]*hours(t)
            Y <- Y + log(1-exp(bCD)*DC0)
            J[2,ib0CD] <- -exp(bCD)*DC0/(1-exp(bCD)*DC0)
            J[2,igCD] <- -exp(bCD)*DC0*hours(t)/(1-exp(bCD)*DC0)
        }
        E <- J %*% LR$cov %*% t(J)
        out[sname,'X'] <- X
        out[sname,'Y'] <- Y
        out[sname,'sX'] <- sqrt(E[1,1])
        out[sname,'sY'] <- sqrt(E[2,2])
        out[sname,'rXY'] <- E[1,2]/(out[sname,'sX']*out[sname,'sY'])
    }
    out
}

plot.calibration <- function(cal,option=1,snames=NULL,i=NULL,...){
    if (is.null(snames)){
        if (stable(cal)) snames <- cal$cal$snames
        else snames <- rownames(cal$cal[[1]]$york)
    }
    if (!is.null(i)) snames <- snames[i]
    if (stable(cal)) calplot_stable(dat=cal,snames=snames,...)
    else calplot_geochronology(dat=cal,option=option,snames=snames,...)
}

calplot_stable <- function(dat,snames=NULL,...){
    num <- dat$m$num
    den <- dat$m$den
    np <- length(num)-1     # number of plot panels
    nr <- ceiling(sqrt(np)) # number of rows
    nc <- ceiling(np/nr)    # number of columns
    oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    for (i in 1:nr){
        for (j in (i+1):max(nc,nr+1)){
            b0names <- names(dat$x)
            xlab <- paste0(num[i],'/',den[i])
            ylab <- paste0(num[j],'/',den[j])
            B <- beta2york(lr=dat,snames=snames,
                           num=num[c(i,j)],den=den[c(i,j)])
            IsoplotR::scatterplot(B,...)
            graphics::mtext(side=1,text=xlab,line=2)
            graphics::mtext(side=2,text=ylab,line=2)
            ell <- IsoplotR::ellipse(dat$cal$x[i],dat$cal$x[j],
                                     dat$cal$cov[c(i,j),c(i,j)])
            graphics::polygon(ell,col='white')
        }
    }
}

calplot_geochronology <- function(dat,option=1,snames=NULL,i=NULL,type,...){
    if (missing(type)) cal <- dat$cal[[1]]
    else cal <- dat$cal[[type]]
    if (is.null(snames)){
        snames <- names(dat$x)
        if (!is.null(i)){
            snames <- snames[i]
        }
    }
    yd <- cal$york[snames,,drop=FALSE]
    num <- cal$num
    den <- cal$den
    xlab <- paste0('log[',num[1],'/',den[1],']')
    ylab <- paste0('log[',num[2],'/',den[2],']')
    if (option==1){
        IsoplotR:::scatterplot(yd,fit=cal$fit,...)
    } else {
        X <- NULL
        Y <- NULL
        for (sname in snames){
            sp <- spot(dat=dat,sname=sname)
            Op <- betapars(spot=sp,ion=num[1])
            Pp <- betapars(spot=sp,ion=den[1])
            Dp <- betapars(spot=sp,ion=num[2])
            Cp <- betapars(spot=sp,ion=num[3])
            newX <- log(Op$sig-Op$bkg) - log(Pp$sig-Pp$bkg) + Op$g*(Op$t-Pp$t)
            newY <- log(Dp$sig-Dp$bkg) - log(Pp$sig-Pp$bkg) + Dp$g*(Dp$t-Pp$t) +
                log(1 - cal$DC0*(Cp$sig-Cp$bkg)/(Dp$sig-Dp$bkg))
            X <- cbind(X,newX)
            Y <- cbind(Y,newY)
        }
        Xlim <- rep(0,2)
        Ylim <- rep(0,2)
        Xlim[1] <- min(yd[snames,'X']-2*yd[snames,'sX'],X)
        Xlim[2] <- max(yd[snames,'X']+2*yd[snames,'sX'],X)
        Ylim[1] <- min(yd[snames,'Y']-2*yd[snames,'sY'],Y)
        Ylim[2] <- max(yd[snames,'Y']+2*yd[snames,'sY'],Y)
        IsoplotR::scatterplot(yd,fit=cal$fit,xlim=Xlim,ylim=Ylim,...)
        matlines(X,Y,lty=1,col='darkgrey')
        if (option>2){
            points(X[1,],Y[1,],pch=21,bg='black')
            points(X[nrow(X),],Y[nrow(Y),],pch=21,bg='white')
        }
    }
    title(caltitle(fit=cal$fit),xlab=xlab,ylab=ylab)
}

caltitle <- function(fit,sigdig=2,type=NA,...){
    args1 <- quote(a%+-%b~'(n='*n*')')
    args2 <- quote(a%+-%b)
    if (is.na(type)){
        intercept <- IsoplotR:::roundit(fit$a[1],fit$a[2],sigdig=sigdig)
        slope <- IsoplotR:::roundit(fit$b[1],fit$b[2],sigdig=sigdig)
        expr1 <- 'slope ='
        expr2 <- 'intercept ='
        list1 <- list(a=slope[1],
                      b=slope[2],
                      u='',
                      n=fit$df+2)
        list2 <- list(a=intercept[1],
                      b=intercept[2],
                      u='')
    }
    call1 <- substitute(e~a,list(e=expr1,a=args1))
    call2 <- substitute(e~a,list(e=expr2,a=args2))
    line1 <- do.call(substitute,list(eval(call1),list1))
    line2 <- do.call(substitute,list(eval(call2),list2))
    line3 <- substitute('MSWD ='~a*', p('*chi^2*')='~b,
                        list(a=signif(fit$mswd,sigdig),
                             b=signif(fit$p.value,sigdig)))
    IsoplotR:::mymtext(line1,line=2,...)
    IsoplotR:::mymtext(line2,line=1,...)
    IsoplotR:::mymtext(line3,line=0,...)
}
