get_alpha <- function(spot){
    nions <- length(spot$ions)
    out <- matrix(0,nions,2)
    colnames(out) <- c('a0','g')
    rownames(out) <- spot$ions
    for (ion in spot$ions){
        tt <- spot$time[,ion]
        sig <- spot$signal[,ion]
        detector <- spot$detector[ion]
        if (spot$nominalblank){
            bkg <- spot$background[detector]
        } else {
            bkg <- spot$signal[,'bkg'] # TODO
        }
        if (spot$type[ion]=='Fc'){
            out[ion,] <- Faraday_blank(tt,sig,bkg)
        } else if (spot$type[detector]=='Em'){
            # TODO add counts to read_data for Em data?
            dwelltime <- spot$dwelltime[ion]
            deadtime <- spot$deadtime[detector]
            edt <- dwelltime - deadtime*sig
            SEM_blank(tt,sig,bkg)
        } else {
            stop('Invalid detector type.')
        }
    }
    out
}

Faraday_blank <- function(tt,sig,bkg){
    a <- log(sig-bkg)
    fit <- lm(a ~ tt)
    out <- fit$coef
    names(out) <- c('a0','g')
    out
}

SEM_blank <- function(tt,sig,bkg){
}
