#' @title read SIMS data
#' @description read a file or folder with ASCII data
#' @param dorf directory or file name
#' @param instrument either \code{'Cameca'} or \code{'SHRIMP'}
#' @param suffix the extension of the data files to be loaded. This
#'     defaults to \code{.asc} if \code{instrument='Cameca'} and
#'     \code{.op} if \code{instrument='SHRIMP'}.
#' @return an object of class \code{simplex}
#' @examples
#' # not run:
#' \dontrun{
#' camdat <- read_data('/path/to/asc/files/',instrument='Cameca',suffix='.asc')
#' plot_timeresolved(camdat[[1]])
#' }
#' @export
read_data <- function(dorf,instrument='Cameca',suffix=NULL){
    if (instrument == 'Cameca') {
        if (is.null(suffix)) suffix <- '.asc'
        out <- read_directory(dorf,instrument='Cameca',suffix=suffix)
    } else if (instrument== 'SHRIMP') {
        out <- read_file(dorf,instrument='SHRIMP',suffix=suffix)
    } else {
        stop('Unsupported instrument')
    }
    class(out) <- c('simplex',instrument)
    out
}

read_directory <- function(dname,instrument='Cameca',suffix='.asc'){
    out <- list()
    fnames <- list.files(dname,pattern=suffix)
    nf <- length(fnames)
    for (i in 1:nf){ # loop through the files
        fname <- fnames[i]
        sname <- tools::file_path_sans_ext(fname)
        out[[sname]] <- read_file(paste0(dname,fname),instrument=instrument)
    }
    out
}

# fname is the complete path to an .asc or .op file
read_file <- function(fname,instrument='Cameca',suffix=NULL){
    suffix <- utils::tail(strsplit(fname,split='[.]')[[1]],n=1)
    if (instrument=='Cameca' & suffix=='asc'){
        out <- read_Cameca_asc(fname)
    } else if (instrument=='SHRIMP' & suffix=='op'){
        out <- read_SHRIMP_op(fname)
    } else if (instrument=='SHRIMP' & suffix=='pd'){
        out <- read_SHRIMP_pd(fname)
    } else {
        stop('Unrecognised file extension.')
    }
    class(out) <- 'spot'
    out
}

read_Cameca_asc <- function(fname){
    ions <- c('Zr90','Zr92','200.5','Zr94',
              'Pb204','Pb206','Pb207','Pb208',
              'U238','ThO2','UO2','270.1')
    f <- file(fname)
    open(f);
    out <- list()
    out$instrument <- 'Cameca'
    while (length(line <- readLines(f,n=1,warn=FALSE)) > 0) {
        if (grepl("ACQUISITION PARAMETERS",line)){
            junk <- readLines(f,n=6,warn=FALSE)
            out$dwelltime <- read_numbers(f,remove=1)
            junk <- readLines(f,n=4,warn=FALSE)
            out$detector <- read_text(f,remove=1)
            names(out$dwelltime) <- ions
            names(out$detector) <- ions
        }
        if (grepl("DETECTOR PARAMETERS",line)) {
            junk <- readLines(f,n=3,warn=FALSE)
            detectors <- NULL
            out$yield <- NULL
            out$background <- NULL
            out$deadtime <- NULL
            while (TRUE){
                line <- readLines(f,n=1,warn=FALSE)
                if (nchar(line)>0){
                    parsed <- parse_line(line)
                    detectors <- c(detectors,parsed[1])
                    out$yield <- c(out$yield, as.numeric(parsed[2]))
                    out$background <- c(out$background, as.numeric(parsed[3]))
                    out$deadtime <- c(out$deadtime, as.numeric(parsed[4]))
                } else {
                    break
                }
            }
            names(out$yield) <- detectors
            names(out$background) <- detectors
            names(out$deadtime) <- detectors
        }
        if (grepl("RAW DATA",line)) {
            cps <- read_asc_block(f,ions=ions)
            out$counts <- round(sweep(cps,MARGIN=2,FUN='*',out$dwelltime))
            out$edt <- effective_dwelltime(out)
        }
        if (grepl("PRIMARY INTENSITY",line)) {
            out$sbm <- read_asc_block(f,ions=ions)
        }
        if (grepl("TIMING",line)) {
            out$time <- read_asc_block(f,ions=ions)
        }
    }
    close(f)
    out
}

read_SHRIMP_op <- function(fname){
    ions <- c('Zr2O','Pb204','bkg','Pb206','Pb207','Pb208','U238','ThO','UO','UO2')
    f <- file(fname)
    open(f);
    out <- list()
    while (TRUE) {
        line <- readLines(f,n=1,warn=FALSE)
        if (length(line)<1){
            break
        } else if (nchar(line)>0){
            spot <- list()
            spot$instrument <- 'SHRIMP'
            sname <- line
            spot$date <- readLines(f,n=1,warn=FALSE)
            spot$set <- read_numbers(f)
            nscans <- read_numbers(f)
            nions <- read_numbers(f)
            spot$deadtime <- 0
            spot$dwelltime <- read_numbers(f)
            names(spot$dwelltime) <- ions
            spot$time <- matrix(0,nscans,nions)
            colnames(spot$time) <- ions
            for (i in 1:nions){
                spot$time[,i] <- read_numbers(f)
            }
            spot$counts <- matrix(0,nscans,nions)
            colnames(spot$counts) <- ions
            for (i in 1:nions){
                spot$counts[,i] <- read_numbers(f)
            }
            spot$edt <- effective_dwelltime(spot)
            spot$sbmbkg <- read_numbers(f)
            spot$sbm <- matrix(0,nscans,nions)
            colnames(spot$sbm) <- ions
            for (i in 1:nions){
                spot$sbm[,i] <- read_numbers(f)
            }
            out[[sname]] <- spot
            junk <- readLines(f,n=1,warn=FALSE)
        }
    }
    close(f)
    out
}

read_SHRIMP_pd <- function(fname){
    ions <- c('Zr2O','Pb204','bkg','Pb206','Pb207','Pb208','U238','ThO','UO','UO2')
    f <- file(fname)
    open(f);
    out <- list()
    while (TRUE) {
        line <- readLines(f,n=1,warn=FALSE)
        if (length(line)<1){
            break
        } else if (nchar(line)>0 & grepl(line,'***',fixed=TRUE)){
            header <- readLines(f,n=4,warn=FALSE)
            spot <- list()
            spot$instrument <- 'SHRIMP'
            namedate <- strsplit(header[[1]],split=', ')[[1]]
            sname <- namedate[1]
            spot$date <- paste(namedate[2:3],collapse=' ')
            spot$set <- split_mixed(header[[2]],1,2)
            nscans <- split_mixed(header[[2]],2,1)
            nions <- split_mixed(header[[2]],3,1)
            spot$sbmbkg <- split_mixed(header[[2]],5,3)
            spot$deadtime <- split_mixed(header[[2]],4,1)
            spot$dwelltime <- read.table(text=readLines(f,n=nions,warn=FALSE))[,4]
            names(spot$dwelltime) <- ions
            spot$time <- matrix(0,nscans,nions)
            spot$counts <- matrix(0,nscans,nions)
            spot$sbm <- matrix(0,nscans,nions)
            block <- readLines(f,n=1+nscans*nions*2,warn=FALSE)[-1]
            colnames(spot$time) <- ions
            colnames(spot$counts) <- ions
            colnames(spot$sbm) <- ions
            for (i in 1:nscans){
                for (j in 1:nions){
                    ii <- (i-1)*nions*2 + (j-1)*2 + 1
                    dat <- read.table(text=block[[ii]])
                    spot$time[i,j] <- as.numeric(dat[3])
                    spot$counts[i,j] <- sum(as.numeric(dat[-(1:4)]))
                    dat <- read.table(text=block[[ii+1]])
                    spot$sbm[i,j] <- sum(as.numeric(dat))
                }
            }
            spot$edt <- effective_dwelltime(spot)
            out[[sname]] <- spot
        }
    }
    close(f)
    out
}

read_text <- function(f,remove=NULL){
    line <- readLines(f,n=1,warn=FALSE)
    parse_line(line,remove=remove)
}
read_numbers <- function(f,remove=NULL){
    parsed <- read_text(f,remove=remove)
    as.numeric(parsed)
}
parse_line <- function(line,remove=NULL){
    parsed <- strsplit(line,'\t|\\s+')[[1]]
    if (is.null(remove)) out <- parsed
    else out <- parsed[-remove]
    out
}
split_mixed <- function(line,i,j,split1=', ',split2=' '){
    chunk <- strsplit(line,split=split1,fixed=TRUE)[[1]][i]
    item <- strsplit(chunk,split=split2)[[1]][j]
    as.numeric(item)
}
read_asc_block <- function(f,ions){
    out <- NULL
    junk <- readLines(f,n=5,warn=FALSE)
    while (TRUE) {
        dat <- read_numbers(f,remove=c(1,2))
        if (length(dat)>0){
            out <- rbind(out,dat)
        } else {
            break
        }
    }
    colnames(out) <- ions
    out
}
# calculate effective dwell time correcting for the dead time
effective_dwelltime <- function(spot){
    if (spot$instrument=='Cameca'){
        deadtime <- spot$deadtime[spot$detector]
    } else if (spot$instrument=='SHRIMP'){
        deadtime <- rep(spot$deadtime,ncol(spot$counts))
    } else {
        stop('Invalid instrument type.')
    }
    lost_time <- sweep(spot$counts,MARGIN=2,FUN='*',deadtime)*1e-9
    -sweep(lost_time,MARGIN=2,FUN='-',spot$dwelltime)
}

subset_samples <- function(dat,prefix='Plesovice',...){
    snames <- names(dat)
    matches <- grep(prefix,snames,...)
    out <- dat[matches]
    class(out) <- class(dat)
    out
}

#' @title define the standards in a dataset
#' @description select a subset of standards from a \code{simplex}
#'     dataset.
#' @param dat an object of class \code{simplex}
#' @param prefix text string to match
#' @param invert logical.  If \code{TRUE} return samples whose names
#'     do _not_ match
#' @param c64 the \eqn{^{206}}Pb/\eqn{^{204}}Pb-ratio of the common Pb
#' @param PbU (optional) true \eqn{^{206}}Pb/\eqn{^{238}}U-ratio of
#'     the age standard
#' @param tst (optional) two-element vector with the age and standard
#'     error of the age standard
#' @return an object of class \code{standard}
#' @examples
#' data(Cameca,package="simplex")
#' stand <- standards(dat=Cameca,prefix='Plesovice')
#' @export
standards <- function(dat,prefix,invert=FALSE,
                      c64=18.7,PbU=NULL,tst=NULL){
    out <- list()
    out$c64 <- c64
    if (is.null(PbU)){
        if (is.null(tst)){
            warning('No standard age or composition was supplied.')
            tst <- Pb76_to_age(dat)
        }
        out$PbU <- IsoplotR:::age_to_Pb206U238_ratio(tt=tst[1],st=tst[2])
    } else {
        out$PbU <- PbU
    }
    out$x <- subset_samples(dat=dat,prefix=prefix,invert=invert)
    class(out) <- 'standard'
    out
}
# get geometric mean Pb207/Pb206 ratio to estimate
# the standard age if not supplied by the user
Pb76_to_age <- function(dat){
    snames <- names(dat)
    ns <- length(snames)
    lPb76 <- rep(0,ns)
    for (i in 1:ns){
        p <- pars(dat[[i]])
        lPb76[i] <- log(sum(p$c7)/sum(p$c6))
    }
    lPb76 <- mean(lPb76)
    slPb76 <- stats::sd(lPb76)/sqrt(ns)
    Pb76 <- exp(lPb76)
    sPb76 <- Pb76*slPb76
    IsoplotR:::get.Pb207Pb206.age.default(x=Pb76,sx=sPb76)
}

#' @title define the samples in a dataset
#' @description select a subset of samples of unknown age from a
#'     \code{simplex} dataset.
#' @param dat an object of class \code{simplex}
#' @param prefix text string to match
#' @param invert logical.  If \code{TRUE} return samples whose names
#'     do _not_ match
#' @return a list of objects of class \code{unknown}
#' @examples
#' data(Cameca,package="simplex")
#' unk <- unknowns(Cameca,prefix='Plesovice',invert=TRUE)
#' @export
unknowns <- function(dat,prefix,invert=FALSE){
    out <- subset_samples(dat=dat,prefix=prefix,invert=invert)
    for (sname in names(out)){
        class(out[[sname]]) <- 'unknown'
    }
    out
}
