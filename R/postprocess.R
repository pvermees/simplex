#' @title convert to a flat data table
#' @description trims the data to an \code{IsoplotR} \code{UPb}
#'     \code{type=5} table
#' @details converts an \code{IsoplotR} \code{type=9} object to a data
#'     table that is formatted as a \code{type=5} input table
#' @param x the output of \code{calibrate}
#' @return a table n \code{IsoplotR} object of class \code{UPb}
#'     (\code{format=5})
#' @examples
#' data(Cameca,package="simplex")
#' stand <- standards(dat=Cameca,prefix='Plesovice')
#' fit <- calibration(stand=stand,oxide='UO2')
#' samp <- unknowns(dat=Cameca,prefix='Qinghu')
#' cal <- calibrate(samp,fit)
#' tab <- data2table(cal)
#' write.csv(tab,file='Qinghu.csv')
#' @export
data2table <- function(x){
    UPb <- logratios2ratios(x)
    ns <- length(x$names)
    nr <- length(x$logratios)
    out <- matrix(0,ns,3*nr)
    rownames(out) <- x$names
    colnames(out) <- c('U238Pb206','s[U238Pb206]',
                       'Pb207Pb206','s[Pb207Pb206]',
                       'Pb204Pb206','s[Pb204Pb206]',
                       'rXY','rXZ','rYZ')
    cormat <- stats::cov2cor(UPb$cov)
    for (i in 1:ns){
        j <- (i-1)*nr+(1:nr)
        out[i,c('U238Pb206','Pb207Pb206',
                'Pb204Pb206')] <- UPb$x[j]
        out[i,c('s[U238Pb206]','s[Pb207Pb206]',
                's[Pb204Pb206]')] <- sqrt(diag(UPb$cov))[j]
        out[i,'rXY'] <- cormat[j[1],j[2]]
        out[i,'rXZ'] <- cormat[j[1],j[3]]
        out[i,'rYZ'] <- cormat[j[2],j[3]]
    }
    out
}
