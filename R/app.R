#setwd('/home/pvermees/Documents/Programming/R/simplex/')
source('R/calibrate.R')
source('R/calibration.R')
source('R/drift.R')
source('R/io.R')
source('R/logratios.R')
source('R/method.R')
source('R/postprocess.R')
source('R/process.R')
source('R/standard.R')
source('R/toolbox.R')
source('R/trim.R')

.simplex <- new.env()

presets <- function(method){
    if (method=='IGG-UPb'){
        load('data/Cameca_UPb.rda')
        out <- Cameca_UPb
    } else if (method=='IGG-UThPb'){
        load('data/Cameca_UThPb.rda')
        out <- Cameca_UThPb
    } else if (method=='IGG-O'){
        load('data/Cameca_oxygen.rda')
        out <- Cameca_oxygen
    } else if (method=='IGG-S'){
        load('data/Cameca_sulphur.rda')
        out <- Cameca_sulphur
    } else if (method=='GA-UPb'){
        load('data/SHRIMP_UPb.rda')
        out <- SHRIMP_UPb
    } else {
        out <- list()
        out$samples <- NULL
        out$method <- defaultmethod(method)
    }
    out
}

as.simplex <- function(dat){
    out <- dat
    ns <- length(dat$samples)
    dc <- "dc"%in%names(out$samples[[1]])
    cameca <- identical(dat$method$instrument,"Cameca")
    for (i in 1:ns){
        if (cameca) {
            names(out$samples[[i]]$detector) <- dat$method$ions
            names(out$samples[[i]]$yield) <- dat$method$detectors
            names(out$samples[[i]]$background) <- dat$method$detectors
        }
        names(out$samples[[i]]$dwelltime) <- dat$method$ions
        names(out$samples[[i]]$dtype) <- dat$method$ions
        names(out$samples[[i]]$deadtime) <- dat$method$detectors
        colnames(out$samples[[i]]$signal) <- dat$method$ions
        colnames(out$samples[[i]]$sbm) <- dat$method$ions
        colnames(out$samples[[i]]$time) <- dat$method$ions
        if (dc) {
            rownames(out$samples[[i]]$dc) <- c('a0','g')
            colnames(out$samples[[i]]$dc) <- dat$method$ions
        }
    }
    class(out) <- "simplex"
    if (dc) class(out) <- append('drift',class(out))
    out
}

driftCorr <- function(x){
    drift(as.simplex(dat=x))
}

driftPlot <- function(x,i){
    plot.drift(x=as.simplex(dat=x),i=as.numeric(i)+1)
}

# f = list of two lists with blocks of text and corresponding filenames
# m = method (currently a string, will be modified to accept lists)
upload <- function(f,m){
    ntcs <- length(f$tcs)
    tcs <- list()
    for (i in 1:ntcs){
        tcs[[f$fns[[i]]]] <- textConnection(f$tcs[[i]])
    }
    read_data(f=tcs,m=m)
}

freeformServer <- function(port=NULL) {
    appDir <- R.utils::getAbsolutePath("inst/www")
    shinylight::slServer(host='0.0.0.0', port=port, appDir=appDir, daemonize=TRUE,
        interface=list(
          presets=presets,
          upload=upload,
          driftCorr=driftCorr,
          driftPlot=driftPlot
        )
    )
}

freeformServer(8000)
