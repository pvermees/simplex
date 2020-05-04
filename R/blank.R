blank <- function(spot,num,den){
    out <- list()
    out$num <- num
    out$den <- den
    nlr <- length(beta$num)
    ntype <- spot$type[beta$num]
    dtype <- spot$type[beta$den]
    Fc <- (ntype%in%'Fc' & dtype%in%'Fc')
    Em <- (ntype%in%'Em' & dtype%in%'Em')
    Mx <- (ntype%in%'Em' & dtype%in%'Fc') | (ntype%in%'Fc' & dtype%in%'Em')
    for (i in 1:nlr){
        numion <- beta$num[i]
        denion <- beta$den[i]
        numtim <- spot$time[,numion]
        dentim <- spot$time[,denion]
        numsig <- spot$sig[,numion]
        densig <- spot$sig[,denion]
        numdet <- spot$detector[numion]
        dendet <- spot$detector[denion]
        numbkg <- get_bkg(spot,numdet)
        denbkg <- get_bkg(spot,dendet)
        if (Fc[i]){
            if (all(numsig>numbkg)){
                numlog <- log(numsig-numbkg)
            } else {
                
            }
            if (all(densig>denbkg)){
                denlog <- log(densig-denbkg)
            } else {
                
            }
            out$lr[,i] <- numlog - denlog
        } else if (Em[i]){
            # TODO
        } else if (Mx[i]){
            # TODO
        } else {
            stop('Invalid detector type.')
        }
    }
    out
}
