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
        simplex <- Cameca_UPb
    } else if (method=='IGG-UThPb'){
        load('data/Cameca_UThPb.rda')
        simplex <- Cameca_UThPb
    } else if (method=='IGG-O'){
        load('data/Cameca_oxygen.rda')
        simplex <- Cameca_oxygen
    } else if (method=='IGG-S'){
        load('data/Cameca_sulphur.rda')
        simplex <- Cameca_sulphur
    } else if (method=='GA-UPb'){
        load('data/SHRIMP_UPb.rda')
        simplex <- SHRIMP_UPb
    } else {
        simplex <- list()
        simplex$samples <- NULL
        simplex$method <- defaultmethod(method)
        class(simplex) <- 'simplex'
    }
    out <- list()
    out$simplex <- simplex
    out$names <- rcnames(simplex)
    out$class <- class(simplex)
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
    out <- restorenames(x$simplex,x$names)
    class(out) <- x$class
    out
}
result2json <- function(x){
    out <- list()
    out$simplex <- x
    out$names <- rcnames(x)
    out$class <- class(x)
    out
}

getdrift <- function(x){
    result2json(drift(x=as.simplex(x)))
}

driftPlot <- function(x,i){
    dat <- as.simplex(x)
    plot.drift(x=dat,i=as.numeric(i)+1)
}

getlogratios <- function(x){
    result2json(logratios(x=as.simplex(x)))
}

logratioPlot <- function(x,i,ratios){
    plot.logratios(x=as.simplex(x),i=as.numeric(i)+1,ratios=ratios)
}

calibrator <- function(x,t=0,option=3,...){
    out <- calibration(as.simplex(x),stand=x$simplex$standard,t=t)
    plot.calibration(out,option=option,...)
    result2json(out)
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
    shinylight::slServer(host='0.0.0.0', port=port,
                         appDir=appDir, daemonize=TRUE,
        interface=list(
          presets=presets,
          upload=upload,
          getdrift=getdrift,
          driftPlot=driftPlot,
          getlogratios=getlogratios,
          logratioPlot=logratioPlot,
          calibrator=calibrator
        )
    )
}

simplex <- function(){
    freeformServer(8000)
}
