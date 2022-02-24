#' @title read or define SIMS acquisition details
#' @description get or set a list of data acquisition properties
#' @param m the name of a data acquisition protocol. Pre-defined
#'     values include \code{'IGG-UPb'}, \code{'GA-UPb'},
#'     \code{'IGG-UThPb'}, \code{'IGG-O'}, and \code{'IGG-S'}).
#' @param instrument one of either \code{SHRIMP} or \code{Cameca}.
#' @param ions vector of labels to be attached to the different ionic
#'     masses that are visited during each sweep.
#' @param num numerator ions of the logratios to be processed in
#'     subsequent data reduction steps.
#' @param den denominator ions of the logratios to be processed in
#'     subsequent data reduction steps.
#' @param blank either (1) a character string with the name of the
#'     channel (e.g. \code{'bkg'}) of the background signal, (2)
#'     \code{'nominal'} a flag indicating that nominal detector
#'     backgrounds recorded in the input file are to be used (only
#'     relevant to Cameca data), or (3) a numerical value representing
#'     a nominal number of background counts (for SHRIMP) or counts
#'     per second (for Cameca).
#' @param description text string with a description of the contents
#' @return an object of class \code{method}
#' @examples
#' fname <- system.file('SHRIMP.op',package='simplex')
#' m <- method(m="GA-UPb",
#'             instrument="SHRIMP",
#'             ions=c('Zr2O','Pb204','bkg','Pb206','Pb207',
#'                    'Pb208','U238','ThO','UO','UO2'),
#'             num=c('Pb204','Pb207','Pb206','UO'),
#'             den=c('Pb206','Pb206','U238','U238'),
#'             blank='bkg',
#'             description="U-Pb dating at Geoscience Australia.")
#' shrimpdat <- read_data(fname,m=m)
#' plot(shrimpdat,i=1)
#' @export
method <- function(m='IGG-UPb',instrument,ions,num,
                   den,blank,description){
    if (m%in%c('IGG-UPb','IGG-UThPb','IGG-O','IGG-S','GA-UPb')){
        out <- defaultmethod(m)
    } else {
        out <- list()
        out$method <- m
    }
    if (!missing(instrument)) out$instrument <- instrument
    if (!missing(ions)) out$ions <- ions
    if (!missing(num)) out$num <- num
    if (!missing(den)) out$den <- den
    if (!missing(blank)) out$bkg <- blank
    if (!missing(description)) out$description <- description
    class(out) <- "method"
    invisible(out)
}

fixmethod <- function(x){
    out <- x
    m <- x$method
    ions <- colnames(x$samples[[1]]$signal)
    if (length(m$ions)!=length(ions)) out$method$ions <- ions
    matches <- all((m$num %in% m$ions) & (m$den %in% m$ions))
    if (!matches) warning("method$num and method$den not found among method$ions")
    invisible(out)
}

defaultmethod <- function(m){
    out <- list()
    out$method <- m
    if (m=='IGG-UPb'){
        out$instrument <- 'Cameca'
        out$ions <- c('Zr90','Zr92','200.5','Zr94',
                      'Pb204','Pb206','Pb207','Pb208',
                      'U238','ThO2','UO2')
        out$num <- c('Pb204','Pb207','Pb206','UO2')
        out$den <- c('Pb206','Pb206','U238','U238')
        out$bkg <- 'nominal'
        out$description <- "Single collector U-Pb dating at CAS-IGG (Beijing)."
    } else if (m=='IGG-UThPb'){
        out$instrument <- 'Cameca'
        out$ions=c('La139','202.5','Pb204','Pb206',
                   'Pb207','Pb208','Th232','U238',
                   'ThO2','UO2')
        out$num <- c('Pb204','Pb207','Pb204','Pb206','UO2','Pb208','ThO2')
        out$den <- c('Pb206','Pb206','Pb208','U238','U238','Th232','Th232')
        out$bkg <- 'nominal'
        out$description <- "Single collector U-Th-Pb dating at CAS-IGG (Beijing)."
    } else if (m=='IGG-O'){
        out$instrument <- 'Cameca'
        out$ions <- c('O16','O17','O18')
        out$num <- 'O18'
        out$den <- 'O16'
        out$bkg <- 'nominal'
        out$description <- "Multicollector oxygen isotope analyses at CAS-IGG (Beijing)."
    } else if (m=='IGG-S'){
        out$instrument='Cameca'
        out$ions <- c('S32','S33','33.96','S34','S36')
        out$num <- c('S33','S34','S36')
        out$den <- c('S32','S32','S32')
        out$bkg <- 'nominal'
        out$description <- "Multicollector sulphur isotope analyses at CAS-IGG (Beijing)."
    } else if (m=='GA-UPb'){
        out$instrument <- 'SHRIMP'
        out$ions <- c('Zr2O','Pb204','bkg','Pb206','Pb207',
                      'Pb208','U238','ThO','UO','UO2')
        out$num <- c('Pb204','Pb207','Pb206','UO')
        out$den <- c('Pb206','Pb206','U238','U238')
        out$bkg <- 'bkg'
        out$description <- "Single collector U-Pb dating at Geoscience Australia."
    } else {
        stop("Invalid method.")
    }
    out
}
