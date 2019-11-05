read_data <- function(dorf,instrument='Cameca',suffix=NULL){
    if (instrument == 'Cameca') {
        if (is.null(suffix)) suffix <- '.asc'
        out <- read_directory(dorf,instrument='Cameca',suffix=suffix)
    } else if (instrument== 'SHRIMP') {
        out <- read_file(dorf,instrument='SHRIMP')
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

# fname is the complete path to an .asc file
read_file <- function(fname,instrument='Cameca'){
    suffix <- tail(strsplit(fname,split='[.]')[[1]],n=1)
    if (instrument=='Cameca' & suffix=='asc'){
        out <- read_Cameca_asc(fname)
    } else if (instrument=='SHRIMP' & suffix=='op'){
        out <- read_SHRIMP_op(fname)
    } else {
        stop('Unrecognised file extension.')
    }
    out
}

read_Cameca_asc <- function(fname){
    ions <- c('Zr90','Zr92','200.5','Zr94',
              'Pb204','Pb206','Pb207','Pb208',
              'U238','ThO2','UO2','270.1')
    f <- file(fname)
    open(f);
    out <- list()
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
        }
        if (grepl("RAW DATA",line)) {
            out$cps <- read_asc_block(f,ions=ions)
            out$counts <- round(sweep(out$cps,MARGIN=2,FUN='*',out$dwelltime))
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
            samp <- list()
            sname <- line
            samp$date <- readLines(f,n=1,warn=FALSE)
            samp$set <- read_numbers(f)
            nscans <- read_numbers(f)
            nions <- read_numbers(f)
            samp$dwelltime <- read_numbers(f)
            samp$time <- matrix(0,nscans,nions)
            colnames(samp$time) <- ions
            for (i in 1:nions){
                samp$time[,i] <- read_numbers(f)
            }
            samp$counts <- matrix(0,nscans,nions)
            for (i in 1:nions){
                samp$counts[,i] <- read_numbers(f)
            }
            samp$sbmbkg <- read_numbers(f)
            samp$sbm <- matrix(0,nscans,nions)
            for (i in 1:nions){
                samp$sbm[,i] <- read_numbers(f)
            }        
            out[[sname]] <- samp
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
    parsed <- strsplit(line,'\t|\\s+')[[1]]
    if (is.null(remove)) out <- parsed
    else out <- parsed[-remove]
    out
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

subset_samples <- function(dat,prefix='Plesovice'){
    snames <- names(dat)
    matches <- grepl(prefix,snames)
    out <- subset(dat,subset=matches)
    class(out) <- class(dat)
    out
}

simplex2isoplotr <- function(alr,format=5){
    snames <- alr$snames
    ns <- length(snames)
    ni <- length(alr$labels)
    ratios <- logratios2ratios(alr)
    out <- matrix(0,ns,9)
    colnames(out) <- c('U238Pb206','s[U238Pb206]',
                       'Pb207Pb206','s[Pb207Pb206]',
                       'Pb204Pb206','s[Pb204Pb206]',
                       'rXY','rXZ','rYZ')
    rownames(out) <- snames
    err <- sqrt(diag(ratios$cov))
    cormat <- stats::cov2cor(ratios$cov)
    for (i in 1:ns){
        i1 <- i      # U8Pb6
        i2 <- i+ns   # Pb76
        i3 <- i+2*ns # Pb46
        out[i,c(1,3,5)] <- ratios$x[c(i1,i2,i3)]
        out[i,c(2,4,6)] <- err[c(i1,i2,i3)]
        out[i,7] <- cormat[i1,i2]
        out[i,8] <- cormat[i1,i3]
        out[i,9] <- cormat[i2,i3]
    }
    out
}
