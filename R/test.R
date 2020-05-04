test_LL_detector <- function(){
    load(file='../data/Cameca_oxygen.rda')
    spot <- Cameca_oxygen[[1]]
    beta <- list()
    beta$num <- c('O17','O18')
    beta$den <- c('O16','O16')
    lr <- log(spot$signal[,beta$num]) - log(spot$signal[,beta$den])
    beta$lr <- matrix(rep(colMeans(lr),times=nrow(spot$signal)),
                      ncol=length(beta$num),byrow=TRUE)
    LL_detector(spot=spot,beta=beta)
}

test_blank <- function(){
    load(file='../data/Cameca_oxygen.rda')
    spot <- Cameca_oxygen[[1]]
    beta <- list()
    beta$num <- c('O17','O18')
    beta$den <- c('O16','O16')
    beta$lr <- log(spot$signal[,beta$num]) - log(spot$signal[,beta$den])
    blank(spot,beta)
}

test_get_a0g <- function(){
    load(file='../data/Cameca_sulfur.rda')
    spot <- Cameca_sulfur[[1]]
    get_a0g(spot)
}

test_alpha <- function(){
    load(file='../data/Cameca_sulfur.rda')
    spot <- Cameca_sulfur[[1]]
    alpha(spot,ions=c('S32','S33','S34','S36'),plot=TRUE)
}
