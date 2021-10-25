setwd('/home/pvermees/Documents/Programming/R/simplex/')
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
    out$names <- rcnames(out$samples)
    out
}

rcnames <- function(dat){
    if (is.list(dat)){
        out <- list()
        for (nm in names(dat)){
            out[[nm]] <- rcnames(dat[[nm]])
        }
    } else if (is.matrix(dat)){
        out <- list(rnames=rownames(dat),cnames=colnames(dat))
    } else {
        out <- names(dat)
    }
    out
}
restorenames <- function(sms,nms){
    out <- sms
    if (is.list(sms)){
        for (nm in names(sms)){
            out[[nm]] <- restorenames(sms[[nm]],nms[[nm]])
        }
    } else if (is.matrix(sms)){
        rownames(out) <- nms$rnames
        colnames(out) <- nms$cnames
    } else {
        names(out) <- nms
    }
    out    
}
as.simplex <- function(x){
    out <- list()
    out$method <- x$method
    out$samples <- restorenames(x$samples,x$names)
    out
}

getdrift <- function(x){
    out <- drift(x=as.simplex(x))
    out$names <- rcnames(out$samples)
    out
}

driftPlot <- function(x,i){
    plot.drift(x=as.simplex(x),i=as.numeric(i)+1)
}

getlogratios <- function(x){
    out <- logratios(as.simplex(x))
    out$names <- rcnames(out$samples)
    out
}

logratioPlot <- function(x,i){
    plot.logratios(x=as.simplex(x),i=as.numeric(i)+1)
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
          getdrift=getdrift,
          driftPlot=driftPlot,
          getlogratios=getlogratios,
          logratioPlot=logratioPlot
        )
    )
}

freeformServer(8000)
