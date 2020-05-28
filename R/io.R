#' @title read SIMS data
#' @description read a file or folder with ASCII data
#' @param dorf directory or file name
#' @param method the name of a data acquisition protocol. One of
#'     \code{instrument='IGG-UPb'}, \code{instrument='GA-zircon'},
#'     \code{instrument='IGG-UThPb'},
#'     \code{instrument='IGG-oxygen'}, or
#'     \code{instrument='IGG-sulfur'}. To create new methods, see
#'     \link{\code{set_method}}.
#' @param suffix the extension of the data files to be loaded. This
#'     defaults to \code{.asc} if \code{instrument='Cameca'} and
#'     \code{.pd} if \code{instrument='SHRIMP'}
#' @return an object of class \code{simplex}
#' @examples
#' # not run:
#' \dontrun{
#' camdat <- read_data('/path/to/asc/files/',instrument='Cameca',suffix='asc')
#' plot_timeresolved(camdat[[1]])
#' }
#' @export
read_data <- function(dorf,suffix,m='IGG-UPb'){
    out <- list()
    if (is.character(m)) m <- method(m)
    if (m$instrument == 'Cameca') {
        if (missing(suffix)) suffix <- 'asc'
        out$samples <- read_directory(dorf,suffix=suffix,m=m)
    } else if (m$instrument== 'SHRIMP') {
        if (missing(suffix)) suffix <- 'pd'
        out$samples <- read_file(dorf,suffix=suffix,m=m)
    } else {
        stop('Unsupported instrument')
    }
    out$method <- m
    class(out) <- 'simplex'
    out
}

read_directory <- function(dname,suffix,m){
    out <- list()
    fnames <- list.files(dname,pattern=suffix)
    nf <- length(fnames)
    for (i in 1:nf){ # loop through the files
        fname <- fnames[i]
        sname <- tools::file_path_sans_ext(fname)
        out[[sname]] <- read_file(paste0(dname,fname),suffix=suffix,m=m)
    }
    out
}

# fname is the complete path to an .asc or .op file
read_file <- function(fname,suffix,m){
    if (m$instrument=='Cameca' & suffix=='asc'){
        out <- read_Cameca_asc(fname=fname,m=m)
    } else if (m$instrument=='SHRIMP' & suffix=='op'){
        out <- read_SHRIMP_op(fname=fname,m=m)
    } else if (m$instrument=='SHRIMP' & suffix=='pd'){
        out <- read_SHRIMP_pd(fname=fname,m=m)
    } else {
        stop('Unrecognised file extension.')
    }
    out
}

read_Cameca_asc <- function(fname,m){
    f <- file(fname)
    open(f);
    out <- list()
    while (length(line <- readLines(f,n=1,warn=FALSE)) > 0) {
        if (grepl("ACQUISITION PARAMETERS",line)){
            junk <- readLines(f,n=6,warn=FALSE)
            out$dwelltime <- read_numbers(f,remove=1)
            junk <- readLines(f,n=4,warn=FALSE)
            out$detector <- read_text(f,remove=1)
            out$dtype <- read_text(f,remove=1)
            names(out$dwelltime) <- m$ions
            names(out$detector) <- m$ions
            names(out$dtype) <- m$ions
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
            out$signal <- read_asc_block(f,ions=m$ions)
        }
        if (grepl("PRIMARY INTENSITY",line)) {
            out$sbm <- read_asc_block(f,ions=m$ions)
        }
        if (grepl("TIMING",line)) {
            out$time <- read_asc_block(f,ions=m$ions)
        }
    }
    close(f)
    out
}

read_SHRIMP_op <- function(fname,m){
    f <- file(fname)
    open(f);
    out <- list()
    while (TRUE) {
        line <- readLines(f,n=1,warn=FALSE)
        if (length(line)<1){
            break
        } else if (nchar(line)>0){
            spot <- list()
            sname <- line
            spot$date <- readLines(f,n=1,warn=FALSE)
            spot$set <- read_numbers(f)
            nscans <- read_numbers(f)
            nions <- read_numbers(f)
            spot$deadtime <- 0
            spot$dwelltime <- read_numbers(f)
            names(spot$dwelltime) <- m$ions
            spot$dtype <- m$dtype
            spot$time <- matrix(0,nscans,nions)
            colnames(spot$time) <- m$ions
            for (i in 1:nions){
                spot$time[,i] <- read_numbers(f)
            }
            spot$signal <- matrix(0,nscans,nions)
            colnames(spot$signal) <- m$ions
            for (i in 1:nions){
                spot$signal[,i] <- read_numbers(f)
            }
            spot$sbmbkg <- read_numbers(f)
            spot$sbm <- matrix(0,nscans,nions)
            colnames(spot$sbm) <- m$ions
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

read_SHRIMP_pd <- function(fname,m){
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
            namedate <- strsplit(header[[1]],split=', ')[[1]]
            sname <- namedate[1]
            spot$date <- paste(namedate[2:3],collapse=' ')
            spot$set <- split_mixed(header[[2]],1,2)
            nscans <- split_mixed(header[[2]],2,1)
            nions <- split_mixed(header[[2]],3,1)
            spot$sbmbkg <- split_mixed(header[[2]],5,3)
            spot$deadtime <- split_mixed(header[[2]],4,1)
            block <- read.table(text=readLines(f,n=nions,warn=FALSE))
            spot$dwelltime <- block[,4]
            names(spot$dwelltime) <- m$ions
            spot$dtype <- rep('Fc',length(m$ions))
            spot$dtype[block[,11]=='COUNTER'] <- 'Em'
            names(spot$dtype) <- m$ions
            spot$time <- matrix(0,nscans,nions)
            spot$signal <- matrix(0,nscans,nions)
            spot$sbm <- matrix(0,nscans,nions)
            block <- readLines(f,n=1+nscans*nions*2,warn=FALSE)[-1]
            colnames(spot$time) <- m$ions
            colnames(spot$signal) <- m$ions
            colnames(spot$sbm) <- m$ions
            for (i in 1:nscans){
                for (j in 1:nions){
                    ii <- (i-1)*nions*2 + (j-1)*2 + 1
                    dat <- read.table(text=block[[ii]])
                    spot$time[i,j] <- as.numeric(dat[3])
                    spot$signal[i,j] <- sum(as.numeric(dat[-(1:4)]))
                    dat <- read.table(text=block[[ii+1]])
                    spot$sbm[i,j] <- sum(as.numeric(dat))
                }
            }
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

subset.simplex <- function(x,prefix=NULL,snames=NULL,i=NULL,...){
    out <- x
    snames <- subset2snames(prefix=prefix,snames=snames,i=i,...)
    out$samples <- x$samples[snames]
    out
}
subset.calibrated <- function(x,prefix=NULL,snames=NULL,i=NULL,...){
    out <- x
    snames <- subset2snames(prefix=prefix,snames=snames,i=i,...)
    out$samples <- x$samples[snames]
    ni <- length(x$calibrated$num)
    i <- which(snames %in% names(x$samples))
    ii <- as.vector(sapply((i-1)*ni,'+',1:ni))
    out$calibrated$lr <- out$calibrated$lr[ii]
    out$calibrated$cov <- out$calibrated$cov[ii,ii]
    out
}
subset2snames <- function(prefix=NULL,snames=NULL,i=NULL,...){
    if (is.null(snames)) snames <- names(dat$samples)
    if (!is.null(i)) snames <- snames[i]
    if (!is.null(prefix)){
        selected <- unlist(lapply(prefix,'grep',snames,...))
        snames <- snames[selected]
    }
    snames
}

spot <- function(dat,sname,i=1,...){
    if (missing(sname)){
        x <- dat$samples[[i]]
        sname <- names(dat)[i]
    } else {
        x <- dat$samples[[sname]]
    }
    out <- dat
    out$samples <- NULL
    out$sname <- sname
    out <- c(out,x)
    class(out) <- 'spot'
    out
}
