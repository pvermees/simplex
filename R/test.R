LL_test <- function(){
    load(file='../data/Cameca_oxygen.rda')
    spot <- Cameca_oxygen[[1]]
    beta <- list()
    beta$num <- c('O17','O18')
    beta$den <- c('O16','O16')
    lr <- log(spot$signal[,beta$num]) - log(spot$signal[,beta$den])
    # single average:
    beta$lr <- matrix(colMeans(lr),ncol=2)
    LL1 <- LL(spot=spot,beta=beta)
    print(paste0('single average: ',LL1))
    # time resolved prediction:
    beta$lr <- matrix(rep(beta$lr,times=nrow(spot$signal)),
                      ncol=length(beta$num),byrow=TRUE)
    LL2 <- LL(spot=spot,beta=beta)
    print(paste0('time resolved: ',LL2))
}
