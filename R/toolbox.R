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

element <- function(ion){
    txt <- gsub("[^a-zA-Z]", "", ion)
    empty <- which(txt %in% "")
    txt[empty] <- ion[empty]
    good <- which(txt %in% elements())
    out <- ion
    out[good] <- txt[good]
    out
}

hours <- function(tt){
    tt/3600
}

background <- function(spot,ions){
    detector <- spot$detector[ions]
    if (spot$nominalblank){
        out <- spot$background[detector]
    } else {
        out <- spot$signal[,'bkg']
    }
    out
}

elements <- function(){
    c('H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg',
      'Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr',
      'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',
      'Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd',
      'Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La',
      'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er',
      'Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au',
      'Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
      'Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md',
      'No','Lr','Rf','Db','Sg','Bh','Hs','Mt')
}

VSMOW <- function(){
    out <- list()
    out$lr <- c(2.00520e-3,0.3799e-3)
    relerr <- c(0.43e-3,1.6e-3)/c(2.00520,0.3799)
    out$cov <- diag(relerr^2)
    labels <- c("O18O16","O17O16")
    names(out$lr) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    out
}

stable <- function(type){
    if (type %in% c("d18O")) return(TRUE)
    else if (type %in% c("U-Pb")) return(FALSE)
    else stop("Invalid data type.")
}

names.simplex <- function(x){
    names(x$x)
}
