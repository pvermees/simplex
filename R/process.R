#' @export
process <- function(f,method,stand,prefix,t=0,exterr=FALSE){
    dat <- read_data(f=f,method=method)
    dc <- drift(x=dat)
    lr <- logratios(x=dc)
    cal <- calibration(lr=lr,stand=stand,prefix=prefix,t=t)
    out <- calibrate(cal,exterr=exterr)
}
