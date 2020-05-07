#' @rdname beta
#' @export
beta <- function(x,...){ UseMethod("beta",x) }
#' @rdname beta
#' @export
beta.default <- function(x,num,den,alpha,...){ stop('No default method.') }
#' @rdname beta
#' @export
beta.spot <- function(x,num,den,alpha,...){
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


