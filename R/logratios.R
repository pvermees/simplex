get_BG <- function(dat,oxide='UO2'){
    snames <- names(dat$x)
    out <- list()
    for (sname in snames){
        spot <- dat$x[[sname]]
        out[[sname]] <- get_bg(spot=spot,oxide=oxide)
    }
    out
}
get_bg <- function(spot,oxide='UO2'){
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

getPbLogRatio <- function(p,cc,num='Pb207',den='Pb206'){
    misfit <- function(b,p,cc,num,den){
        bij <- b + log(p[[num]]$d) - log(p[[den]]$d)
        LL <- p[[num]]$n*bij - (p[[num]]$n+p[[den]]$n)*log(1+exp(bij))
        -sum(LL)
    }
    init <- mean(cc$bdc76)
    optimise(misfit,interval=init+c(-5,5),
             p=p,cc=cc,num=num,den=den)$minimum
}
