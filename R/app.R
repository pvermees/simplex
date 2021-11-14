presets <- function(method){
    if (method=='IGG-UPb'){
        data("Cameca_UPb",package="simplex")
        simplex <- Cameca_UPb
    } else if (method=='IGG-UThPb'){
        data("Cameca_UThPb",package="simplex")
        simplex <- Cameca_UThPb
    } else if (method=='IGG-O'){
        data("Cameca_oxygen",package="simplex")
        simplex <- Cameca_oxygen
    } else if (method=='IGG-S'){
        data("Cameca_sulphur",package="simplex")
        simplex <- Cameca_sulphur
    } else if (method=='GA-UPb'){
        data("SHRIMP_UPb",package="simplex")
        simplex <- SHRIMP_UPb
    } else {
        simplex <- list()
        simplex$samples <- NULL
        simplex$method <- defaultmethod(method)
        class(simplex) <- 'simplex'
    }
    result2json(simplex)
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
    out$multi <- multicollector(x)
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
    out <- as.simplex(x)
    plot.drift(x=out,i=as.numeric(i)+1)
    result2json(out)
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
    dat <- as.simplex(x)
    if (x$fixedslope) slope <- x$slope
    else slope <- NULL
    out <- calibration(dat,stand=x$simplex$standard,slope=slope)
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
    cal <- calibrate_it(x)
    if (identical(x$IsoplotRformat,'U-Pb')){
        method <- 'U-Pb'
    } else if (identical(x$IsoplotRformat,'U-Th-Pb')){
        method <- 'U-Th-Pb'
    } else if (identical(x$IsoplotRformat,'Th-Pb')){
        method <- 'Th-Pb'
    }
    out <- simplex2IsoplotR(cal,method=method)
    as.list(as.data.frame(out$x))
}

# f = list of two lists with blocks of text and corresponding filenames
upload <- function(f,x){
    ntcs <- length(f$tcs)
    tcs <- list()
    for (i in 1:ntcs){
        tcs[[f$fns[[i]]]] <- textConnection(f$tcs[[i]])
    }
    result2json(read_data(f=tcs,m=as.simplex(x)$method))
}

freeformServer <- function(port=NULL,host='127.0.0.1',test=FALSE) {
    if (test) appDir <- R.utils::getAbsolutePath("inst/www")
    else appDir <- system.file("www",package="simplex")
    shinylight::slServer(host=host,port=port,appDir=appDir,
                         daemonize=!is.null(port),
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
#' @param port The port on which to listen. If not provided, a random
#'     unused port will be chosen and a browser window opened on that
#'     port.
#' @examples
#' \donttest{simplex()}
#' @export
simplex <- function(port=NULL){
    freeformServer(port)
}
