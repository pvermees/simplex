#' @title read SIMS data
#' @description read ASCII data files
#' @param f file name(s), may include wildcards (\code{*.asc},
#'     \code{*.op} or \code{*.pd}).
#' @param m an object of class \code{method} OR the name of a data
#'     acquisition protocol (one of \code{'IGG-UPb'}, \code{'GA-UPb'},
#'     \code{'IGG-UThPb'}, \code{'IGG-O'}, or \code{'IGG-S'}). To
#'     create new methods, see \code{method}.
#' @return an object of class \code{simplex}
#' @examples
#' fname <- system.file('SHRIMP.op',package='simplex')
#' shrimpdat <- read_data(fname,m='GA-UPb')
#' plot(shrimpdat,i=1)
#' @export
read_data <- function(f,m='IGG-UPb'){
    out <- list()
    if (is.character(m)){
        out$method <- method(m=m)
    } else {
        out$method <- m
    }
    if ("textConnection" %in% class(f[[1]])){
        out$samples <- read_samples_tc(tc=f,m=out$method)
    } else {
        out$samples <- read_samples_fn(fn=f,m=out$method)
    }
    class(out) <- 'simplex'
    fixmethod(out)
}
read_samples_tc <- function(tc,m){
    out <- list()
    ntc <- length(tc)
    fn <- names(tc)
    for (i in 1:ntc){
        if (m$instrument == 'Cameca') {
            sname <- tools::file_path_sans_ext(fn[i])
            out[[sname]] <- read_file(tc[[i]],m=m)
        } else if (m$instrument== 'SHRIMP') {
            ext <- tools::file_ext(fn[i])
            out <- c(out,read_file(tc[[i]],m=m,ext=ext))
        } else {
            stop('Unsupported instrument')
        }
    }
    out
}
read_samples_fn <- function(fn,m){
    out <- list()
    for (fname in Sys.glob(fn)){
        f <- file(fname)
        open(f);
        if (m$instrument == 'Cameca') {
            sname <- tools::file_path_sans_ext(fname)
            out[[sname]] <- read_file(f,m=m)
        } else if (m$instrument== 'SHRIMP') {
            ext <- tools::file_ext(fname)
            out <- c(out,read_file(f,m=m,ext=ext))
        } else {
            stop('Unsupported instrument')
        }
        close(f)
    }
    out
}

read_file <- function(f,m,ext=NA){
    if (m$instrument=='Cameca'){
        out <- read_Cameca_asc(f=f,mions=m$ions)
    } else if (m$instrument=='SHRIMP') {
        if (identical(ext,'op')){
            out <- read_SHRIMP_op(f=f,mions=m$ions)
        } else if (identical(ext,'pd')) {
            out <- read_SHRIMP_pd(f=f,mions=m$ions)
        } else {
            stop('Invalid file extension')
        }
    } else {
        stop('Unrecognised file extension.')
    }
    out
}

read_Cameca_asc <- function(f,mions=NULL){
    out <- list()
    while (length(line <- readLines(f,n=1,warn=FALSE)) > 0) {
        if (grepl("X POSITION",line)){
            parsed <- parse_line(line)
            out$x <- as.numeric(parsed[2])
            out$y <- as.numeric(parsed[4])
        }
        if (grepl("ACQUISITION PARAMETERS",line)){
            junk <- readLines(f,n=1,warn=FALSE)
            ions <- read_text(f,remove=1)
            if (length(ions)==length(mions)) ions <- mions
            junk <- readLines(f,n=4,warn=FALSE)
            out$dwelltime <- read_numbers(f,remove=1)
            junk <- readLines(f,n=4,warn=FALSE)
            out$detector <- read_text(f,remove=1)
            out$dtype <- read_text(f,remove=1)
            names(out$dwelltime) <- ions
            names(out$detector) <- ions
            names(out$dtype) <- ions
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
            out$signal <- read_asc_block(f,ions=ions)
        }
        if (grepl("PRIMARY INTENSITY",line)) {
            out$sbm <- read_asc_block(f,ions=ions)
        }
        if (grepl("TIMING",line)) {
            out$time <- read_asc_block(f,ions=ions)
        }
    }
    out
}

read_SHRIMP_op <- function(f,mions=NULL){
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
            if (nions==length(mions)) ions <- mions
            else ions <- paste0('m',1:nions)
            names(spot$dwelltime) <- ions
            spot$dtype <- rep('Em',length(ions))
            names(spot$dtype) <- ions
            spot$time <- matrix(0,nscans,nions)
            colnames(spot$time) <- ions
            for (i in 1:nions){
                spot$time[,i] <- read_numbers(f)
            }
            spot$signal <- matrix(0,nscans,nions)
            colnames(spot$signal) <- ions
            for (i in 1:nions){
                spot$signal[,i] <- read_numbers(f)
            }
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
    out
}

read_SHRIMP_pd <- function(f,mions=NULL){
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
            block <- utils::read.table(text=readLines(f,n=nions,warn=FALSE))
            ions <- block[,1]
            if (length(ions)==length(mions)) ions <- mions
            spot$dwelltime <- block[,4]
            names(spot$dwelltime) <- ions
            spot$detector <- block[,11]
            spot$dtype <- rep('Fc',length(ions))
            spot$dtype[spot$detector=='COUNTER'] <- 'Em'
            names(spot$dtype) <- ions
            spot$time <- matrix(0,nscans,nions)
            spot$signal <- matrix(0,nscans,nions)
            spot$sbm <- matrix(0,nscans,nions)
            block <- readLines(f,n=1+nscans*nions*2,warn=FALSE)[-1]
            colnames(spot$time) <- ions
            colnames(spot$signal) <- ions
            colnames(spot$sbm) <- ions
            for (i in 1:nscans){
                for (j in 1:nions){
                    ii <- (i-1)*nions*2 + (j-1)*2 + 1
                    dat <- utils::read.table(text=block[[ii]])
                    spot$time[i,j] <- as.numeric(dat[3])
                    spot$signal[i,j] <- sum(as.numeric(dat[-(1:4)]))
                    dat <- utils::read.table(text=block[[ii+1]])
                    spot$sbm[i,j] <- sum(as.numeric(dat))
                }
            }
            out[[sname]] <- spot
        }
    }
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
    parsed <- trimws(strsplit(line,'\t')[[1]])
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
    snames <- subset2snames(dat=x,prefix=prefix,snames=snames,i=i,...)
    out$samples <- x$samples[snames]
    out
}
subset.calibrated <- function(x,prefix=NULL,snames=NULL,i=NULL,...){
    out <- x
    snames <- subset2snames(dat=x,prefix=prefix,snames=snames,i=i,...)
    out$samples <- x$samples[snames]
    ni <- length(x$calibrated$num)
    i <- which(names(x$samples) %in% snames)
    ii <- as.vector(sapply((i-1)*ni,'+',1:ni))
    out$calibrated$lr <- out$calibrated$lr[ii]
    out$calibrated$cov <- out$calibrated$cov[ii,ii]
    out
}
subset2snames <- function(dat,prefix=NULL,snames=NULL,i=NULL,...){
    if (is.null(snames)) snames <- names(dat$samples)
    if (!is.null(i)) snames <- snames[i]
    if (!is.null(prefix)){
        selected <- unlist(lapply(prefix,'grep',snames,...))
        snames <- snames[selected]
    }
    snames
}

spot <- function(dat,sname=NULL,i=1,...){
    if (is.null(sname)){
        x <- dat$samples[[i]]
        sname <- names(dat$samples)[i]
    } else {
        x <- dat$samples[[sname]]
    }
    out <- dat
    out$samples <- NULL
    out$sname <- sname
    if (!is.null(out$outliers)) out$outliers <- dat$outliers[[sname]]
    out <- c(out,x)
    class(out) <- 'spot'
    out
}

#' @title plot simplex data
#' @description plot time resolved SIMS data
#' @param x an object of class \code{simplex}
#' @param sname the sample name to be shown
#' @param i the sample number to be shown
#' @param ... optional arguments to be passed on to the generic
#'     \code{plot} function.
#' @examples
#' data('SHRIMP_UPb',package='simplex')
#' plot(SHRIMP_UPb,i=1)
#' @method plot simplex
#' @export
plot.simplex <- function(x,sname=NULL,i=1,...){
    plot.spot(x=spot(dat=x,sname=sname,i=i),...)
}
#' @export
plot.spot <- function(x,...){
    ions <- names(x$dwelltime)
    np <- length(ions)      # number of plot panels
    nr <- ceiling(sqrt(np)) # number of rows
    nc <- ceiling(np/nr)    # number of columns
    oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    for (ion in ions){
        graphics::plot(x$time[,ion],x$signal[,ion],
                       type='p',xlab='',ylab='',...)
        graphics::mtext(side=1,text='t',line=2)
        graphics::mtext(side=2,text=ion,line=2)
    }
    graphics::par(oldpar)
}
