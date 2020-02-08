get_bg <- function(dat,oxide='UO2'){
    snames <- names(dat$x)
    ns <- length(snames)
    init <- matrix(0,ns,2)
    colnames(init) <- c('b','g')
    rownames(init) <- snames
    out <- list(O=init,U238=init,blank=init,
                Pb204=init,Pb206=init,Pb207=init)
    for (sname in snames){
        spot <- dat$x[[sname]]
        bg <- get_bg_spot(spot=spot,oxide=oxide)
        out$O[sname,] <- bg['O',]
        out$U238[sname,] <- bg['U238',]
        out$Pb206[sname,] <- bg['Pb206',]
        out$Pb207[sname,] <- bg['Pb207',]
        out$Pb204[sname,] <- bg['Pb204',]
        out$blank[sname,] <- bg['blank',]
    }
    out
}
get_bg_spot <- function(spot,oxide='UO2'){
    p <- pars(spot,oxide=oxide)
    out <- matrix(0,6,2)
    colnames(out) <- c('b','g')
    rownames(out) <- c('O','U238','Pb204','Pb206','Pb207','blank')
    out['O',] <- bg_helper(p=p$O)
    out['U238',] <- bg_helper(p=p$U238)
    out['Pb206',] <- bg_helper(p=p$Pb206)
    out['Pb207',] <- b_helper(p=p$Pb207,g=out['Pb206','g'])
    out['Pb204',] <- b_helper(p=p$Pb204,g=out['Pb206','g'])
    if (sum(p$blank$n)>0)
        out['blank','b'] <- log(sum(p$blank$n)) - log(sum(p$blank$d))
    else
        out['blank','b'] <- -Inf
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

b_helper <- function(p,g){
    if (sum(p$n)>0)
        b <- log(sum(p$n)) - log(sum(exp(g*p$t)*p$d))
    else
        b <- -Inf
    c(b,g)
}

blank_correct <- function(bg,sname,tt,mass){
    b <- bg[[mass]][sname,'b']
    g <- bg[[mass]][sname,'g']
    bb <- bg$blank[sname,'b']
    db <- bb-b-g*tt
    if (all(db<0)) out <- log(1-exp(db))
    else out <- -Inf
    out
}

drift_correct <- function(num,den,g){
    log(num$c/den$c) + g*(den$t-num$t)
}
