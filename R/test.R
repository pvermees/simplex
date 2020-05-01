test_LL_detector <- function(){
    load(file='../data/Cameca_oxygen.rda')
    spot <- Cameca_oxygen[[1]]
    beta <- list()
    beta$num <- c('O17','O18')
    beta$den <- c('O16','O16')
    lr <- log(spot$signal[,beta$num]) - log(spot$signal[,beta$den])
    # single average:
    beta$lr <- matrix(rep(colMeans(lr),times=nrow(spot$signal)),
                      ncol=length(beta$num),byrow=TRUE)
    LL_detector(spot=spot,beta=beta)
}

test_LL_blank <- function(){
    load(file='../data/Cameca_oxygen.rda')
    spot <- Cameca_oxygen[[1]]
    get_alpha(spot)
}
