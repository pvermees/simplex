get_bg <- function(dat,oxide='UO2'){
    snames <- names(dat$x)
    ns <- length(snames)
    out <- matrix(0,ns,6)
    colnames(out) <- c('bO','bU','b6','gO','gU','gPb')
    rownames(out) <- snames
    for (sname in snames){
        spot <- dat$x[[sname]]
        p <- pars(spot,oxide=oxide)
        out[sname,c('bO','gO')] <- bg_helper(tt=p$tO,dd=p$dO,nn=p$nO)
        out[sname,c('bU','gU')] <- bg_helper(tt=p$tU,dd=p$dU,nn=p$nU)
        out[sname,c('b6','gPb')] <- bg_helper(tt=p$t6,dd=p$d6,nn=p$n6)
    }
    out
}

bg_helper <- function(tt,nn,dd){
    misfit <- function(bg,tt,nn,dd){
        LL <- nn*(bg[1]+bg[2]*tt+log(dd)) - exp(bg[1]+bg[2]*tt)*dd
        -sum(LL)
    }
    init <- c(0,0)
    fit <- optim(init,misfit,tt=tt,nn=nn,dd=dd)
    fit$par
}
