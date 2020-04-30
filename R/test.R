LL_test <- function(){
    load(file='../data/Cameca.rda')
    spot <- Cameca[[1]]
    beta <- list()
    LL(spot=spot,beta=beta)
}
