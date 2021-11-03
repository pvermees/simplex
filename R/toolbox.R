# calculate effective dwell time correcting for the dead time
effective_dwelltime <- function(spot){
    if (spot$method$instrument=='Cameca'){
        deadtime <- spot$deadtime[spot$detector]
    } else if (spot$method$instrument=='SHRIMP'){
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

isotope <- function(ion){
    as.numeric(gsub("[^0-9]", "", ion))
}

hours <- function(tt){
    tt/3600
}

seconds <- function(tt){
    tt*3600
}

background <- function(spot,ions){
    detector <- spot$detector[ions]
    if (spot$method$nominalblank){
        out <- spot$background[detector]
    } else if (spot$m$instrument=='Cameca'){
        out <- spot$signal[,'bkg']
    } else if (spot$m$instrument=='SHRIMP'){
        out <- spot$signal[,'bkg']/spot$dwelltime['bkg']
    } else {
        stop('Illegal background option')
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
    out$lr <- log(c(0.3799e-3,2.00520e-3))
    relerr <- c(1.6e-3,0.43e-3)/c(0.3799,2.00520)
    out$cov <- diag(relerr^2)
    labels <- c("O17O16","O18O16")
    names(out$lr) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    out
}

troilite <- function(){
    out <- list()
    S2346 <- c(126.948,22.6436,6515)
    out$lr <- -log(S2346)
    relerr <- c(0.047,0.0020,20)/S2346
    out$cov <- diag(relerr^2)
    labels <- c("S33S32","S34S32","S36S32")
    names(out$lr) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    out    
}

stable <- function(dat){
    type <- datatype(dat)
    if (type %in% c("oxygen")) return(TRUE)
    else if (type %in% c("sulphur")) return(TRUE)
    else if (type %in% c("U-Pb")) return(FALSE)
    else if (type %in% c("U-Th-Pb")) return(FALSE)
    else stop("Invalid data type.")
}

datatype <- function(x){
    ions <- x$method$ions
    if (all(c("U238","Th232","Pb204",
              "Pb206","Pb207","Pb208")%in%ions) &&
               any(c("UO","UO2")%in%ions) &&
               any(c("ThO","ThO2")%in%ions)){
        out <- "U-Th-Pb"
    } else if (all(c("U238","Pb206","Pb206","Pb207")%in%ions) &&
        any(c("UO","UO2")%in%ions)){
        out <- "U-Pb"
    } else if (all(c("O16","O17","O18")%in%ions)){
        out <- "oxygen"
    } else if (all(c('S32','S33','S34','S36')%in%ions)){
        out <- "sulphur"
    }
    out
}

chronometer <- function(x,type){
    dt <- datatype(x)
    if (identical(dt,"U-Pb"))
        out <- "U238-Pb206"
    else if (identical(dt,"U-Th-Pb") & type %in% c(1,'UPb'))
        out <- "U238-Pb206"
    else if (identical(dt,"U-Th-Pb") & type %in% c(2,'ThPb'))
        out <- "Th232-Pb208"
    else
        stop('Invalid chronometer')
    out
}

faraday <- function(spot,ion=NULL){
    if (is.null(ion)) out <- any(spot$dtype[spot$method$ions]=='Fc')
    else out <- spot$dtype[ion]=='Fc'
    out
}
