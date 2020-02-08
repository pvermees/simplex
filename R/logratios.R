get_bg <- function(dat,oxide='UO2'){
    snames <- names(dat$x)
    ns <- length(snames)
    init <- matrix(0,ns,2)
    colnames(init) <- c('b','g')
    rownames(init) <- snames
    out <- list(O=init,U238=init,Pb206=init,blank=init)
    for (sname in snames){
        spot <- dat$x[[sname]]
        p <- pars(spot,oxide=oxide)
        out$O[sname,] <- bg_helper(p=p$O)
        out$U238[sname,] <- bg_helper(p=p$U238)
        out$Pb206[sname,] <- bg_helper(p=p$Pb206)
        if (sum(p$blank$n)>0)
            out$blank[sname,'b'] <- log(sum(p$blank$n)) - log(sum(p$blank$d))
        else
            out$blank[sname,'b'] <- -Inf
    }
    out
}

bg_helper <- function(p){
    misfit <- function(par,p){
        b <- par[1]
        g <- par[2]
        LL <- p$n*(b+g*p$t+log(p$d)) - exp(b+g*p$t)*p$d
        -sum(LL)
    }
    init <- c(0,0)
    fit <- optim(par=init,fn=misfit,p=p)
    fit$par
}

blank_correction <- function(bg,bb,tt){
    log(1-exp(bb-bg['b']-bg['g']*tt))
}
