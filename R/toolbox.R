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
        out <- as.numeric(blk)
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

# Optimise with some fixed parameters 
# Like optim, but with option to fix some parameters.
# parms: Parameters to potentially optimize in fn
# fixed: A vector of TRUE/FALSE values indicating which parameters in
#   parms to hold constant (not optimize). If TRUE, the corresponding
#   parameter in fn() is fixed. Otherwise it's variable and optimised over.
# lower, upper: a vector with the lower and upper ends of the search
#   ranges of the free parameters, respectively
# Originally written by Barry Rowlingson, modified by PV
optifix <- function(parms, fixed, fn, gr = NULL, ...,
                    method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"), 
                    lower = -Inf, upper = Inf, control = list(), hessian = FALSE){
    force(fn)
    force(fixed) 
    .npar <- length(parms)
    .pnames <- names(parms)
    .fixValues <- parms[fixed]
    names(.fixValues) <- names(parms)[fixed]
    .parStart <- parms[!fixed]
    names(.parStart) <- names(parms)[!fixed]
  
    .fn <- function(parms,pnames=names(parms),...){
        .par <- rep(NA,.npar)
        .par[!fixed] <- parms
        .par[fixed] <- .fixValues
        names(.par) <- .pnames
        fn(.par,...)
    }

    if(!is.null(gr)){
        .gr <- function(parms,pnames=names(parms),...){
            .gpar <- rep(NA,.npar)
            .gpar[!fixed] <- parms
            .gpar[fixed] <- .fixValues
            names(.gpar) <- .pnames
            gr(.gpar,...)[!fixed]
        }
    } else {
        .gr <- NULL
    }

    .opt <- stats::optim(.parStart,.fn,.gr,...,method=method,
                        lower=lower,upper=upper,
                        control=control,hessian=hessian) 
    
    .opt$fullpars <- rep(NA,.npar) 
    .opt$fullpars[fixed] <- .fixValues 
    .opt$fullpars[!fixed] <- .opt$par
    names(.opt$fullpars) <- names(parms)
    .opt$fixed <- fixed

    # remove fullpars (PV)
    .opt$par <- .opt$fullpars
    .opt$fullpars <- NULL
    
    return(.opt)
}
