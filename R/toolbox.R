# calculate effective dwell time correcting for the dead time
effective_dwelltime <- function(spot){
    if (spot$instrument=='Cameca'){
        deadtime <- spot$deadtime[spot$detector]
    } else if (spot$instrument=='SHRIMP'){
        deadtime <- rep(spot$deadtime,ncol(spot$signal))
    } else {
        stop('Invalid instrument type.')
    }
    lost_time <- sweep(spot$signal,MARGIN=2,FUN='*',deadtime)*1e-9
    -sweep(lost_time,MARGIN=2,FUN='-',spot$dwelltime)
}
