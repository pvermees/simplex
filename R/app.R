presets <- function(method){
    if (method=='IGG-UPb'){
        simplex <- get(load('data/Cameca_UPb.rda'))
    } else if (method=='IGG-UThPb'){
        simplex <- get(load('data/Cameca_UThPb.rda'))
    } else if (method=='IGG-O'){
        simplex <- get(load('data/Cameca_oxygen.rda'))
    } else if (method=='IGG-S'){
        simplex <- get(load('data/Cameca_sulphur.rda'))
    } else if (method=='GA-UPb'){
        simplex <- get(load('data/SHRIMP_UPb.rda'))
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

getdatatype <- function(x){
    dat <- as.simplex(x)
    datatype(dat)
}

getdrift <- function(x){
    out <- drift_helper(x=as.simplex(x),gui=TRUE)
    result2json(out)
}

driftPlot <- function(x,i){
    dat <- as.simplex(x)
    plot.drift(x=dat,i=as.numeric(i)+1)
}

getlogratios <- function(x){
    out <- logratios_helper(x=as.simplex(x),gui=TRUE)
    result2json(out)
}

logratioPlot <- function(x,i,ratios){
    plot.logratios(x=as.simplex(x),i=as.numeric(i)+1,ratios=ratios)
}

getstandard <- function(preset){
    standard(preset)
}

calibrator <- function(x,...){
    out <- calibration(as.simplex(x),stand=x$simplex$standard)
    plot.calibration(out,...)
    result2json(out)
}

calibrate_it <- function(x){
    dat <- as.simplex(x)
    selection <- subset(dat,prefix=x$sampleprefix)
    out <- calibrate(selection)
}

calibrateSamples <- function(x){
    out <- calibrate_it(x)
    plot.calibrated(out)
    out
}

plotresults <- function(x){
    cal <- calibrate_it(x)
    if (stable(cal)){
        d <- delta(cal)
        plot.delta(d)
    } else {
        UPb <- simplex2IsoplotR(cal)
        IsoplotR::concordia(UPb)
    }
}

resultstable <- function(x){
    cal <- calibrate_it(x)
    if (stable(cal)){
        d <- delta(cal)
    } else {
        d <- cal
    }
    tab <- data2table(d)
    rownames(tab) <- NULL
    as.data.frame(tab)
}

export2isoplotr <- function(x){
    if (identical(x$IsoplotRformat,'U-Pb')){
        method <- 'U-Pb'
    } else if (identical(x$IsoplotRformat,'U-Th-Pb')){
        method <- 'U-Th-Pb'
    } else if (identical(x$IsoplotRformat,'Th-Pb')){
        method <- 'Th-Pb'
    }
    out <- simplex2IsoplotR(as.simplex(x),method=method)
    as.list(as.data.frame(out$x))
}

# f = list of two lists with blocks of text and corresponding filenames
# m = method (currently a string, will be modified to accept lists)
upload <- function(f,m){
    ntcs <- length(f$tcs)
    tcs <- list()
    for (i in 1:ntcs){
        tcs[[f$fns[[i]]]] <- textConnection(f$tcs[[i]])
    }
    result2json(read_data(f=tcs,m=m))
    
}

freeformServer <- function(port=NULL) {
    appDir <- R.utils::getAbsolutePath("inst/www")
    shinylight::slServer(host='0.0.0.0', port=port,
                         appDir=appDir, daemonize=TRUE,
        interface=list(
            presets=presets,
            upload=upload,
            getdatatype=getdatatype,
            getdrift=getdrift,
            driftPlot=driftPlot,
            getlogratios=getlogratios,
            logratioPlot=logratioPlot,
            getstandard=getstandard,
            calibrator=calibrator,
            calibrateSamples=calibrateSamples,
            plotresults=plotresults,
            resultstable=resultstable,
            export2isoplotr=export2isoplotr
        )
    )
}

#' @title Graphical User Interface for simplex
#' @description Shinylight app for simplex
#' @examples
#' \donttest{simplex()}
#' @export
simplex <- function(){
    freeformServer()
}
