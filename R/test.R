LL_test <- function(){
    load(file='../data/Cameca_oxygen.rda')
    spot <- Cameca_oxygen[[1]]
    beta <- list()
    beta$num <- c('O17','O18')
    beta$den <- c('O16','O16')
    lr <- colMeans(log(spot$signal[,beta$num]) -
                   log(spot$signal[,beta$den]))
    beta$lr <- matrix(lr,ncol=2)
    LL(spot=spot,beta=beta)
}
