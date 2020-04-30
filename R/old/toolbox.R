pars <- function(spot,parent='U238',daughter='Pb206',oxide='UO2'){
    init <- function(spot,mass){
        out <- list()
        out$t <- hours(spot$time[,mass])
        out$n <- spot$counts[,mass]
        out$d <- spot$edt[,mass]
        out$c <-  out$n/out$d
        out
    }
    CamecaBlank <- function(spot,mass){
        out <- list()
        bkg <- spot$background[spot$detector]
        names(bkg) <- names(spot$detector)
        out$t <- spot$time[,mass]
        nt <- length(out$t)
        out$c <- rep(bkg[mass],nt)
        out$d <- rep(1,nt)
        out$n <- out$c*out$d
        out
    }
    SHRIMPblank <- function(spot){
        out <- list()
        out$n <- spot$counts[,'bkg']
        out$t <- hours(spot$time[,'bkg'])
        out$d <- spot$edt[,'bkg']
        out$c <- out$n/out$d
        out
    }
    out <- list(
        Pb204 = init(spot,'Pb204'),
        Pb206 = init(spot,'Pb206'),
        Pb207 = init(spot,'Pb207'),
        Pb208 = init(spot,'Pb208'),
        oxide = oxide,
        parent = parent,
        daughter = daughter
    )
    out[[parent]] <- init(spot,parent)
    out[[oxide]] <- init(spot,oxide)
    if (spot$instrument=='Cameca'){
        out$blank <- CamecaBlank(spot,mass=parent)
    } else {
        out$blank <- SHRIMPblank(spot)
    }
    out
}

hours <- function(tt){
    tt/3600
}
