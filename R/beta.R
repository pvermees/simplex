#' @rdname beta
#' @export
beta <- function(x,...){ UseMethod("beta",x) }
#' @rdname beta
#' @export
beta.default <- function(x,...){ stop('No default method.') }
#' @rdname beta
#' @export
beta.spot <- function(x,num,den,a,...){
    out <- matrix(0,2,length(num))
    rownames(out) <- c('b0','g')
    colnames(out) <- paste(num,den,sep='/')
    groups <- groupbypairs(num,den)
    for (gr in groups){
        ni <- length(gr$num)
        if (gr$elements['num']==gr$elements['den']){
            init <- rep(0,ni)
            fit <- optim(par=init,f=SS_b0,method='BFGS',spot=x,
                         num=gr$num,den=gr$den,a=a)
        } else {
            init <- rep(0,ni+1)
            fit <- optim(par=init,f=SS_b0g,method='BFGS',spot=x,
                         num=gr$num,den=gr$den,a=a)
        }
    }
    out
}

SS_b0g <- function(b0g,spot,num,den,a){
    ni <- length(num)
    nt <- nrow(spot$time)
    out <- 0
    for (i in 1:ni){
        a0 <- rep(0,nt) # TODO
        fit <- optim(a0,SS_a0,b0g=b0g,spot=spot,
                     num=num[i],den=den[i],a=a)
        out <- out + fit$value
    }
    out
}
SS_b0 <- function(b0,spot,num,den,a){
    SS_b0g(c(b0,0),spot=spot,num=num,den=den,a=a)
}

SS_a0 <- function(a0,b0g,spot,num,den,a){
    nb0 <- length(b0g)-1
    b0 <- b0g[1:nb0]
    g <- b0g[nb0+1]
    Ni <- spot$signal[,num]
    Nb <- background(spot,num)
    Ng <- a['g',num]
    Nt <- spot$time[,num]
    Di <- spot$signal[,den]
    Db <- background(spot,den)
    Dg <- a['g',den]
    Dt <- spot$time[,den]
    SS <- (Ni-Nb-exp(a0+b0+g*Dt+Ng*(Nt-Dt)))^2 + (Di-Db-exp(a0))^2
    sum(SS)
}

groupbypairs <- function(num,den){
    out <- list()
    pairs <- elementpairs(num,den)
    nele <- element(num)
    dele <- element(den)
    for (i in 1:ncol(pairs)){
        matches <- which((nele %in% pairs[1,i]) & (dele %in% pairs[2,i]))
        out[[i]] <- list()
        out[[i]]$elements <- pairs[,i]
        names(out[[i]]$elements) <- c('num','den')
        out[[i]]$num <- num[matches]
        out[[i]]$den <- den[matches]
    }
    out
}

elementpairs <- function(num,den){
    nele <- element(num)
    dele <- element(den)
    out <- rbind(nele[1],dele[1])
    if (length(num)>1){
        for (i in 2:length(num)){
            alreadyin <- (nele[i] %in% out[1,] || dele[i] %in% out[2,])
            if (!alreadyin){
                out <- cbind(out,c(nele[i],dele[i]))
            }
        }
    }
    out
}
