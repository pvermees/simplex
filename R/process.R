#' @title process SIMS data
#' @description one stop shop for SIMS data reduction, including
#'     including the drift correction, logratio calculation, and
#'     calibration.
#' @param f file name(s), may include wildcards (\code{*.asc},
#'     \code{*.op} or \code{*.pd}).
#' @param method an object of class \code{method} OR the name of a
#'     data acquisition protocol (one of \code{'IGG-UPb'},
#'     \code{'GA-UPb'}, \code{'IGG-UThPb'}, \code{'IGG-O'}, or
#'     \code{'IGG-S'}). To create new methods, see \code{method}.
#' @param stand an object of class \code{standard}
#' @param t the analysis time to which the logratio signals should be
#'     regressed.
#' @param exterr include the uncertainty associated with the
#'     standard calibration in the error propagation?
#' @return an object of class \code{calibrated}
#' @examples
#' \dontrun{
#' m <- method('GA-UPb')
#' s <- standard(preset="Temora",prefix='TEM')
#' fname <- system.file('SHRIMP.pd',package='simplex')
#' cd <- process(f=fname,m=m,stand=s)
#' tab <- data2table(cd)
#' }
#' @export
process <- function(f,m,stand,t=0,exterr=FALSE){
    dat <- read_data(f=f,m=m)
    dc <- drift(x=dat)
    lr <- logratios(x=dc)
    cal <- calibration(lr=lr,stand=stand,t=t)
    calibrate(cal,exterr=exterr)
}
