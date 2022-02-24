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
    blk <- spot$method$bkg
    numeric <- suppressWarnings(!is.na(as.numeric(blk)))
    if (identical(blk,'nominal')){
        out <- spot$background[detector]
    } else if (numeric){
        if (spot$m$instrument=='Cameca'){
            out <- as.numeric(blk)
        } else if (spot$m$instrument=='SHRIMP'){
            out <- as.numeric(blk)/spot$dwelltime[ions]
        } else {
            stop('Illegal background option')
        }
    } else {
        if (spot$m$instrument=='Cameca'){
            out <- spot$signal[,blk]
        } else if (spot$m$instrument=='SHRIMP'){
            out <- spot$signal[,blk]/spot$dwelltime[blk]
        } else {
            stop('Illegal background option')
        }
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

multicollector <- function(x,...){ UseMethod("multicollector",x) }
multicollector.default <- function(x,...){
    detectors <- x$detector
    numion <- length(detectors)
    numdet <- length(unique(detectors))
    if (numdet==numion) return(TRUE)
    else if (numdet==1) return(FALSE)
    else stop('Simplex does not currently manage mixed ',
              'peak-hopping and multidetection runs.')    
}
multicollector.simplex <- function(x){
    multicollector.default(x$samples[[1]])
}
multicollector.spot <- function(x){
    multicollector.default(x)
}

datatype <- function(x){
    ions <- x$method$ions
    if (all(c("U238","Pb206","Pb206","Pb207")%in%ions) &&
        any(c("UO","UO2")%in%ions)){
        out <- "U-Pb"
    } else if (all(c("Th232","Pb204","Pb208")%in%ions) &&
               any(c("ThO","ThO2")%in%ions)){
        out <- "Th-Pb"
    } else if (all(c("O16","O17","O18")%in%ions)){
        out <- "oxygen"
    } else if (all(c('S32','S33','S34','S36')%in%ions)){
        out <- "sulphur"
    }
    out
}

chronometer <- function(x){
    dt <- datatype(x)
    if (identical(dt,"U-Pb"))
        out <- "U238-Pb206"
    else if (identical(dt,"Th-Pb"))
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

outlier <- function(x,sname=NULL,i=1,j){
    if (is.null(sname)) sname <- names(x$samples)[i]
    if (missing(j)){
        out <- x$outliers[[sname]]
    } else {
        out <- x
        if (is.null(x$outliers)){
            # does not yet have any outliers -> create
            out$outliers <- list()
            out$outliers[[sname]] <- j
        } else if (is.null(x$outliers[[sname]])){
            # no outliers for this sample -> add
            out$outliers[[sname]] <- j
        } else if (j %in% x$outliers[[sname]]){
            # does not yet have this outlier -> add
            out$outliers[[sname]] <- c(out$outliers[[sname]],j)
        } else {
            # already has this outlier -> remove
            jj <- which(x$outliers[[sname]] %in% j)
            out$outliers[[sname]] <- out$outliers[[sname]][-jj]
        }
    }
    out
}
