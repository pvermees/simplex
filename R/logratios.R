get_B0G <- function(dat,parent='U238',oxide='UO2'){
    snames <- names(dat)
    out <- list()
    for (sname in snames){
        spot <- dat[[sname]]
        out[[sname]] <- get_b0g(spot=spot,parent=parent,oxide=oxide)
    }
    out
}
get_b0g <- function(spot,parent='U238',oxide='UO2'){
    p <- pars(spot,parent=parent,oxide=oxide)
    out <- matrix(0,7,2)
    colnames(out) <- c('b0','g')
    rownames(out) <- c(oxide,parent,'Pb204','Pb206','Pb207','Pb208','blank')
    out[oxide,] <- b0g_helper(p=p[[p$oxide]])
    out[parent,] <- b0g_helper(p=p[[p$parent]])
    out['Pb206',] <- b0g_helper(p=p$Pb206)
    out['Pb207',] <- b0_helper(p=p$Pb207,g=out['Pb206','g'])
    out['Pb208',] <- b0_helper(p=p$Pb208,g=out['Pb206','g'])
    out['Pb204',] <- b0_helper(p=p$Pb204,g=out['Pb206','g'])
    if (sum(p$blank$n)>0)
        out['blank','b0'] <- log(sum(p$blank$n)) - log(sum(p$blank$d))
    else
        out['blank','b0'] <- -Inf
    out
}

b0g_helper <- function(p){
    misfit <- function(par,p){
        b0 <- par[1]
        g <- par[2]
        LL <- p$n*(b0+g*p$t+log(p$d)) - exp(b0+g*p$t)*p$d
        -sum(LL)
    }
    init <- c(0,0)
    fit <- optim(par=init,fn=misfit,p=p)
    fit$par
}

b0_helper <- function(p,g){
    if (sum(p$n)>0)
        b0 <- log(sum(p$n)) - log(sum(exp(g*p$t)*p$d))
    else
        b0 <- -Inf
    c(b0,g)
}

getPbLogRatios <- function(dat){
    snames <- names(dat)
    ns <- length(snames)
    out <- list()
    out$snames <- snames
    out$num <- c('Pb204','Pb207','Pb208')
    out$den <- c('Pb206','Pb206','Pb206')
    out$x <- rep(0,3*ns)
    out$cov <- matrix(0,3*ns,3*ns)
    for (i in 1:ns){
        spot <- dat[[i]]
        p <- pars(spot=spot)
        b0g <- get_b0g(spot=spot)
        lr <- getPbLogRatio(p,b0g)
        cormat <- cov2cor(lr$cov)
        j <- c(0,ns,2*ns)+i
        out$x[j] <- lr$x
        out$cov[j,j] <- sqrt(diag(lr$cov))
    }
    out
}
getPbLogRatio <- function(p,b0g){
    misfit_helper <- function(b,p,b0g,num,den){
        bpc <- b + A2Corr(p=p,b0g=b0g,num=num,den=den)
        bpn <- bpc + log(p[[num]]$d) - log(p[[den]]$d)
        LL <- LLbinom(bn=bpn,nnum=p[[num]]$n,nden=p[[den]]$n)
        sum(LL)
    }
    misfit <- function(par,p,b0g){
        LL46 <- misfit_helper(b=par[1],p=p,b0g=b0g,num='Pb204',den='Pb206')
        LL76 <- misfit_helper(b=par[2],p=p,b0g=b0g,num='Pb207',den='Pb206')
        LL86 <- misfit_helper(b=par[3],p=p,b0g=b0g,num='Pb208',den='Pb206')
        -(LL46 + LL76 + LL86)
    }
    init46 <- log(mean(p$Pb204$c)) - log(mean(p$Pb206$c))
    init76 <- log(mean(p$Pb207$c)) - log(mean(p$Pb206$c))
    init86 <- log(mean(p$Pb208$c)) - log(mean(p$Pb206$c))
    init <- c(init46,init76,init86)
    fit <- stats::optim(par=init,fn=misfit,p=p,b0g=b0g,hessian=TRUE)
    out <- list()
    out$x <- fit$par
    out$cov <- solve(fit$hessian)
    out
}

# atomic to cps correction: 
# to be added to atomic and subtracted from cps logratios
A2Corr <- function(p,b0g,num='Pb207',den='Pb206'){
    tnum <- p[[num]]$t
    tden <- p[[den]]$t
    # 1. drift correction
    if (identical(den,'U238')){
        dc <- b0g[num,'g']*(tnum-tden)
    } else {
        dc <- b0g[num,'g']*tnum - b0g[den,'g']*tden
    }
    # 2. blank correction
    bcnum <- blank_correct(b0g=b0g,tt=tnum,mass=num)
    bcden <- blank_correct(b0g=b0g,tt=tden,mass=den)
    dc - bcnum + bcden
}

blank_correct <- function(b0g,tt,mass){
    b0 <- b0g[mass,'b0']
    g <- b0g[mass,'g']
    bb <- b0g['blank','b0']
    db <- bb-b0-g*tt
    if (all(db<0)) out <- log(1-exp(db))
    else out <- rep(-Inf,length(tt))
    out
}

# bn = logratio of counts, nnum = counts of num, nden = counts of den
LLbinom <- function(bn,nnum,nden){
    nnum*bn - (nnum+nden)*log(1+exp(bn))
}
