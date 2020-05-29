process <- function(dorf,meth,stand,prefix,t=0,exterr=FALSE){
    dat <- read_data(dorf,m=meth)
    dc <- drift(x=dat)
    lr <- logratios(x=dc)
    cal <- calibration(lr=lr,stand=stand,prefix=prefix,t=t)
    out <- calibrate(cal,exterr=exterr)
}
