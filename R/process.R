process <- function(dof,meth,stand,prefix,suffix,t=0,exterr=FALSE){
    dat <- read_data(dof,suffix=suffix,m=meth)
    dc <- drift(x=dat)
    lr <- logratios(x=dc)
    cal <- calibration(lr=lr,stand=stand,prefix=prefix,t=t)
    out <- calibrate(cal,exterr=exterr)
}
