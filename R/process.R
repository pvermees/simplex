process <- function(dof,method){
    dat <- read_data(datadir,method=method,suffix=suffix)
    dc <- drift(x=dat,ions=ions)
    lr <- logratios(x=dc,num=num,den=den)
    cal <- calibration(lr=lr,stand=st,prefix=stand)
    cd <- calibrate(cal,exterr=TRUE)
}
