presets <- function(method){
    if (method=='IGG-UPb'){
        simplex <- utils::data('Cameca_UPb',package='simplex',
                               envir=environment())
    } else if (method=='IGG-UThPb'){
        simplex <- utils::data('Cameca_UThPb',package='simplex',
                               envir=environment())
    } else if (method=='IGG-O'){
        simplex <- utils::data('Cameca_oxygen',package='simplex',
                               envir=environment())
    } else if (method=='IGG-S'){
        simplex <- utils::data('Cameca_sulphur',package='simplex',
                               envir=environment())
    } else if (method=='GA-UPb'){
        simplex <- utils::data('SHRIMP_UPb',package='simplex',
                               envir=environment())
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

preset2standard <- function(x){
    preset <- x$calibration$preset
    measured <- identical(x$calibration$standtype,'measured')
    out <- as.simplex(x)
    out$calibration$stand <- standard(preset=preset,measured=measured)
    result2json(out)
}

t2stand <- function(x){
    tst <- x$calibration$tst
    measured <- identical(x$calibration$standtype,'measured')
    out <- as.simplex(x)
    out$calibration$stand <- standard(tst=tst,measured=measured)
    result2json(out)
}
d2stand <- function(x){
    d <- x$calibration$del
    delval <- as.numeric(d$delval[1,])
    refval <- as.numeric(d$refval[1,])
    nr <- length(delval)
    delcov <- matrix(0,nr,nr)
    refcov <- matrix(0,nr,nr)
    for (r in 1:nr){
        delcov[r,] <- as.numeric(d$delcov[[r]])
        refcov[r,] <- as.numeric(d$refcov[[r]])
    }
    names(delval) <- d$ratios
    rownames(delcov) <- d$ratios
    colnames(delcov) <- d$ratios
    names(refval) <- d$ratios
    rownames(refcov) <- d$ratios
    colnames(refcov) <- d$ratios
    del <- list(val=delval,cov=delcov)
    ref <- list(val=refval,cov=refcov)
    out <- as.simplex(x)
    out$calibration$stand <- standard(del=del,ref=ref)
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
    if (identical(x$calibration$caltype,"average")) prng <- NULL
    else prng <- dat$calibration$pairing
    prfx <- dat$calibration$prefix
    if (length(x$standards)>0){
        snms <- x$standards
    } else {
        snms <- NULL
    }
    out <- calibration(dat,stand=stnd,pairing=prng,prefix=prfx,snames=snms)
    plot.calibration(out,show.numbers=x$shownum,...)
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
    plot.calibrated(out,show.numbers=x$shownum)
    result2json(out)
}

calibratedTable <- function(x){
    cal <- calibrate_it(x)
    tab <- data2table.calibrated(cal,log=x$log)
    rownames(tab) <- NULL
    as.data.frame(tab)
}

preset2deltaref <- function(x){
    ref <- x$delta$preset
    if (identical(ref,'other')){
        dat <- as.simplex(x)
        out <- dat$calibration$stand$ref
    } else {
        out <- do.call(ref,list())
    }
    out$ratios <- names(out$val)
    out
}

convert2delta <- function(x){
    # 1. re-calibrate
    dat <- calibrate_it(x)
    # 2. record reference
    val <- x$delta$val
    nv <- length(val)
    names(val) <- x$delta$ratios
    ref <- list(preset=x$delta$preset,val=val)
    # 3. get delta values
    del <- delta(dat,ref=ref,log=identical(x$delta$type,'delta-prime'))
    tab <- data2table.delta(del)
    rownames(tab) <- NULL
    as.data.frame(tab)
}

convert2IsoplotR <- function(x){
    dat <- calibrate_it(x)
    isodat <- simplex2IsoplotR(dat,method=x$IsoplotRtype)
    as.data.frame(isodat$x)
}

export2IsoplotR <- function(x){
    as.list(convert2IsoplotR(x))
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
            preset2standard=preset2standard,
            createcalibration=createcalibration,
            createpairing=createpairing,
            t2stand=t2stand,
            d2stand=d2stand,
            calibrator=calibrator,
            calibrateSamples=calibrateSamples,
            calibratedTable=calibratedTable,
            preset2deltaref=preset2deltaref,
            convert2delta=convert2delta,
            convert2IsoplotR=convert2IsoplotR,
            export2IsoplotR=export2IsoplotR
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
