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
data2table <- function(lr){
    r <- logratios2ratios(lr)
    ns <- length(lr$snames)
    nc <- length(lr$x)/ns
    out <- matrix(0,ns,nc)
    for (i in 1:nc){
    }
    r
}

logratios2ratios <- function(lr){
    out <- list()
    out$snames <- lr$snames
    out$x <- exp(-lr$x)
    J <- -diag(out$x)
    out$cov <- J %*% lr$cov %*% t(J)
    out
}
