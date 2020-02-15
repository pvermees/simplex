get_B0G <- function(dat,oxide='UO2'){
    snames <- names(dat$x)
    out <- list()
    for (sname in snames){
        spot <- dat$x[[sname]]
        out[[sname]] <- get_b0g(spot=spot,oxide=oxide)
    }
    out
}
get_b0g <- function(spot,oxide='UO2'){
    p <- pars(spot,oxide=oxide)
    out <- matrix(0,6,2)
    colnames(out) <- c('b0','g')
    rownames(out) <- c('O','U238','Pb204','Pb206','Pb207','blank')
    out['O',] <- b0g_helper(p=p$O)
    out['U238',] <- b0g_helper(p=p$U238)
    out['Pb206',] <- b0g_helper(p=p$Pb206)
    out['Pb207',] <- b0_helper(p=p$Pb207,g=out['Pb206','g'])
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

# atomic to cps correction: 
# to be added to atomic and subtracted from cps logratios
A2Corr <- function(p,b0g,num='Pb207',den='Pb206'){
    dcnum <- b0g[num,'g']*p[[num]]$t # drift correction for num
    dcden <- b0g[den,'g']*p[[den]]$t # drift correction for den
    bcnum <- blank_correct(b0g=b0g,tt=p[[num]]$t,mass=num)
    bcden <- blank_correct(b0g=b0g,tt=p[[den]]$t,mass=den)
    dcnum - dcden - bcnum + bcden
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
