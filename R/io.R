read_data <- function(dorf,instrument='Cameca',suffix=NULL){
    if (instrument == 'Cameca') {
        if (is.null(suffix)) suffix <- '.asc'
        out <- read_directory(dorf,instrument='Cameca',suffix=suffix)
    } else if (instrument== 'SHRIMP') {
        out <- read_file(dorf,instrument='SHRIMP')
    } else {
        stop('Unsupported instrument')
    }
}

read_directory <- function(dname,instrument='Cameca',suffix='.asc'){
    out <- list()
    class(out) <- c('simplex',instrument)
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
    f <- file(fname)
    open(f);
    out <- list()
    while (length(line <- readLines(f, n=1, warn=FALSE)) > 0) {
        if (grepl("ACQUISITION PARAMETERS",line)) {
            block <- readLines(f, n=12, warn=FALSE)
            ions <- unlist(strsplit(block[2], split='\t'))[-1]
            out$ions <- gsub(' ','',ions)
            out$dwelltime <- as.numeric(unlist(strsplit(block[7], split='\t'))[-1])
            detector <- unlist(strsplit(block[12], split='\t'))[-1]
            out$detector <- gsub(' ','',detector)
            names(out$dwelltime) <- out$ions
            names(out$detector) <- out$ions
        }
        if (grepl("DETECTOR PARAMETERS",line)) {
            block <- readLines(f, n=3, warn=FALSE)
            detectors <- NULL
            out$yield <- NULL
            out$background <- NULL
            while (TRUE){
                line <- readLines(f, n=1, warn=FALSE)
                if (grepl("CORRECTION FACTORS",line)){
                    break
                } else if (nchar(line)>0){
                    detectorpars <- gsub(' ','',unlist(strsplit(line, split='\t')))
                    detectors <- c(detectors,detectorpars[1])
                    out$yield <- c(out$yield, as.numeric(detectorpars[2]))
                    out$background <- c(out$background, as.numeric(detectorpars[3]))
                }
            }
            names(out$yield) <- detectors
            names(out$background) <- detectors
        }
        if (grepl("RAW DATA",line)) {
            out$cps <- NULL
            junk <- readLines(f, n = 5, warn = FALSE)
            while ((line <- readLines(f, n=1, warn=FALSE)) != '') {
                dat <- as.numeric(unlist(strsplit(line, split='\t')))[-c(1,2)]
                out$cps <- rbind(out$cps,dat)
            }
            colnames(out$cps) <- out$ions
            out$counts <- round(sweep(out$cps,MARGIN=2,FUN='*',out$dwelltime))
        }
        if (grepl("PRIMARY INTENSITY",line)) {
            out$sbm <- NULL
            junk <- readLines(f, n = 5, warn = FALSE)
            while ((line <- readLines(f, n=1, warn=FALSE)) != '') {
                dat <- as.numeric(unlist(strsplit(line, split='\t')))[-c(1,2)]
                out$sbm <- rbind(out$sbm,dat)
            }
            colnames(out$sbm) <- out$ions
        }
        if (grepl("TIMING",line)) {
            out$time <- NULL
            junk <- readLines(f, n = 5, warn = FALSE)
            while (length(line <- readLines(f, n=1, warn=FALSE)) > 0) {
                dat <- as.numeric(unlist(strsplit(line, split='\t')))[-c(1,2)]
                out$time <- rbind(out$time,dat)
            }
            colnames(out$time) <- out$ions
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
        print(line)
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

read_numbers <- function(f){
    line <- readLines(f,n=1,warn=FALSE)
    parsed <- strsplit(line,'\t| ')[[1]]
    as.numeric(parsed)
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
