#' @title read or define SIMS acquisition details
#' @description get or set a list of data acquisition properties
#' @param method the name of a data acquisition protocol. Pre-defined
#'     values include \code{'IGG-UPb'}, \code{'GA-UPb'},
#'     \code{'IGG-UThPb'}, \code{'IGG-O'}, and \code{'IGG-S'}).
#' @param instrument one of either \code{SHRIMP} or \code{Cameca}.
#' @param ions vector of labels to be attached to the different ionic
#'     masses that are visited during each sweep.
#' @param num numerator ions of the logratios to be processed in
#'     subsequent data reduction steps.
#' @param den denominator ions of the logratios to be processed in
#'     subsequent data reduction steps.
#' @param oxide the ion label to be used as a uranium or thorium oxide
#'     reference for U-Pb and U-Th-Pb calibration.
#' @param nominalblank logical flag indicating whether the blank is
#'     measured during every sweep, or whether to use a nominal blank
#'     value for each detector.
#' @param description text string with a description of the contents
#' @return an object of class \code{method}
#' @examples
#' fname <- system.file('SHRIMP.pd',package='simplex')
#' shrimpdat <- read_data(fname,method='GA-UPb')
#' plot(shrimpdat,i=1)
#' @export
method <- function(method='IGG-UPb',instrument,ions,
                   num,den,oxide,nominalblank,description){
    if (method%in%c('IGG-UPb','IGG-UThPb','IGG-O','IGG-S','GA-UPb')){
        out <- defaultmethod(method)
    } else {
        out <- list()
        out$method <- method
    }
    if (!missing(instrument)) out$instrument <- instrument
    if (!missing(ions)) out$ions <- ions
    if (!missing(num)) out$num <- num
    if (!missing(den)) out$den <- den
    if (!missing(oxide)) out$oxide <- oxide
    if (!missing(nominalblank)) out$nominalblank <- nominalblank
    if (!missing(description)) out$description <- description
    class(out) <- "method"
    invisible(out)
}

defaultmethod <- function(m){
    out <- list()
    out$method <- m
    if (m=='IGG-UPb'){
        out$instrument <- 'Cameca'
        out$ions <- c('Zr90','Zr92','200.5','Zr94',
                      'Pb204','Pb206','Pb207','Pb208',
                      'U238','ThO2','UO2','270.1')
        out$num <- c('Pb204','Pb207','Pb208','Pb206','UO2')
        out$den <- c('Pb206','Pb206','Pb206','U238','U238')
        out$oxide <- c(U='UO2')
        out$nominalblank <- TRUE
        out$description <- "Single collector U-Pb dating at CAS-IGG (Beijing)."
    } else if (m=='IGG-UThPb'){
        out$instrument='Cameca'
        out$ions=c('La139','202.5','Pb204','Pb206',
                   'Pb207','Pb208','Th232','U238',
                   'ThO2','UO2')
        out$num <- c('Pb204','Pb207','Pb208','Pb206','UO2','Pb208','ThO2')
        out$den <- c('Pb206','Pb206','Pb206','U238','U238','Th232','Th232')
        out$oxide <- c(U='UO2',Th='ThO2')
        out$nominalblank <- TRUE
        out$description <- "Single collector U-Th-Pb dating at CAS-IGG (Beijing)."
    } else if (m=='IGG-O'){
        out$instrument <- 'Cameca'
        out$ions <- c('O16','O17','O18')
        out$num <- c('O17','O18')
        out$den <- c('O16','O16')
        out$nominalblank <- TRUE
        out$description <- "Multicollector oxygen isotope analyses at CAS-IGG (Beijing)."
    } else if (m=='IGG-S'){
        out$instrument='Cameca'
        out$ions <- c('S32','S33','33.96','S34','S36')
        out$num <- c('S33','S34','S36')
        out$den <- c('S32','S32','S32')
        out$nominalblank <- TRUE
        out$description <- "Multicollector sulphur isotope analyses at CAS-IGG (Beijing)."
    } else if (m=='GA-UPb'){
        out$instrument <- 'SHRIMP'
        out$ions <- c('Zr2O','Pb204','bkg','Pb206','Pb207',
                      'Pb208','U238','ThO','UO','UO2')
        out$dtype <- rep('Em',length(out$ions))
        out$num <- c('Pb204','Pb207','Pb208','Pb206','UO')
        out$den <- c('Pb206','Pb206','Pb206','U238','U238')
        out$oxide <- c(U='UO')
        out$nominalblank <- FALSE
        out$description <- "Single collector U-Pb dating at Geoscience Australia."
    } else {
        stop("Invalid method.")
    }
    out
}
