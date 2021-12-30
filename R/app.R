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

driftPlot <- function(x){
    dat <- as.simplex(x)
    i <- as.numeric(x$i)+1
    out <- drift(dat,i=i)
    plot.drift(x=out,i=i)
    result2json(out)
}

getlogratios <- function(x){
    out <- logratios_helper(x=as.simplex(x),gui=TRUE)
    result2json(out)
}

logratioPlot <- function(x,ratios){
    dat <- as.simplex(x)
    i <- as.numeric(x$i)+1
    dc <- drift(dat,i=i)
    out <- logratios(dat,i=i)
    plot.logratios(x=out,i=i,ratios=ratios)
    result2json(out)
}

logratioTable <- function(x){
    tab <- data2table.logratios(as.simplex(x),log=x$log,addxy=x$xy)
    rownames(tab) <- NULL
    as.data.frame(tab)
}

getstandard <- function(x){
    preset <- x$calibration$preset
    measured <- identical(x$calibration$caltype,'measured')
    out <- as.simplex(x)
    out$calibration$stand <- standard(preset,measured)
    result2json(out)
}

createcalibration <- function(x){
    measured <- (x$calibration$standtype=='measured')
    dat <- as.simplex(x)
    s <- skeletonstand(dat,measured=measured)
    p <- pairing(dat,s)
    out <- dat
    out$calibration <- list(stand=s,pairing=p)
    result2json(out)
}
createpairing <- function(x){
    dat <- as.simplex(x)
    s <- dat$calibration$stand
    p <- pairing(dat,s)
    out <- dat
    out$calibration <- list(stand=s,pairing=p)
    result2json(out)
}

calibrator <- function(x,...){
    dat <- as.simplex(x)
    stnd <- dat$calibration$stand
    prng <- dat$calibration$pairing
    prfx <- dat$calibration$prefix
    if (length(x$standards)>0){
        snms <- x$standards
    } else {
        snms <- NULL
    }
    out <- calibration(dat,stand=stnd,pairing=prng,prefix=prfx,snames=snms)
    plot.calibration(out,...)
    result2json(out)
}

calibrate_it <- function(x){
    dat <- as.simplex(x)
    if (length(x$samples)>0){
        snms <- x$samples
    } else {
        snms <- NULL
    }
    selection <- subset(dat,prefix=x$sampleprefix,snames=snms)
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
upload <- function(f,m){
    ntcs <- length(f$tcs)
    tcs <- list()
    for (i in 1:ntcs){
        tcs[[f$fns[[i]]]] <- textConnection(f$tcs[[i]])
    }
    out <- read_data(f=tcs,m=m)
    result2json(out)
}

freeformServer <- function(port=NULL,host='127.0.0.1',
                           test=FALSE,daemonize=!is.null(port)) {
    if (test) appDir <- R.utils::getAbsolutePath("inst/www")
    else appDir <- system.file("www",package="simplex")
    shinylight::slServer(host=host,port=port,appDir=appDir,daemonize=daemonize,
        interface=list(
            presets=presets,
            upload=upload,
            getdatatype=getdatatype,
            getdrift=getdrift,
            driftPlot=driftPlot,
            getlogratios=getlogratios,
            logratioPlot=logratioPlot,
            logratioTable=logratioTable,
            getstandard=getstandard,
            createcalibration=createcalibration,
            createpairing=createpairing,
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
