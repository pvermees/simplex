beta <- function(spot,num,den){
    out <- list()
    out$num <- num
    out$den <- den
    numele <- element(num)
    denele <- element(den)
    nlr <- length(num)
    for (i in 1:nlr){
        numion <- num[i]
        denion <- den[i]
        numtim <- spot$time[,numion]
        dentim <- spot$time[,denion]
        numsig <- spot$sig[,numion]
        densig <- spot$sig[,denion]
        numdet <- spot$detector[numion]
        dendet <- spot$detector[denion]
        numbkg <- get_bkg(spot,numdet)
        denbkg <- get_bkg(spot,dendet)
    }
    
}

get_a0g <- function(spot){
    el <- element(spot$ions)
    EL <- unique(el)
    out <- matrix(0,nrow=length(spot$ions),ncol=2)
    colnames(out) <- c('a0','g')
    for (e in EL){
        i <- which(el %in% e)
        ions <- spot$ions[i]
        out[i,] <- get_a0g_helper(spot,ions)
    }
    out
}

# groups all the ions from the same element
get_a0g_helper <- function(spot,ions){
    init_g <- function(tt,sb){
        imax <- which.max(colSums(sb))
        lsb <- log(sb[,imax])
        lm(lsb ~ tt)$coef[2]
    }
    init_g <- function(tt,sig,bkg){
        
    }
    tt <- spot$time[,ions,drop=FALSE]
    sig <- spot$signal[,ions,drop=FALSE]
    bkg <- get_bkg(spot=spot,ions=ions)
    g <- init_g(tt=tt,sb=sig-bkg)
    a <- init_a(tt=tt,sb=sig,bkg=bkg)
}

