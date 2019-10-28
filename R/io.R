# dname is a directory containing .asc files
read_directory <- function(dname,instrument='Cameca',suffix=NULL){
    if (instrument == 'Cameca') {
        if (is.null(suffix)) suffix <- '.asc'
    } else if (instrument== 'SHRIMP') {
        if (is.null(suffix)) suffix <- '.op'
    } else {
        stop('Unsupported instrument')
    }
    out <- list()
    class(out) <- instrument
    fnames <- list.files(dname,pattern=suffix)
    nf <- length(fnames)
    for (i in 1:nf){ # loop through the files
        fname <- fnames[i]
        sname <- tools::file_path_sans_ext(fname)
        out[[sname]] <- read_file(paste0(dname,fname),
                                  instrument=instrument,suffix=suffix)
    }
    out
}

# fname is the complete path to an .asc file
read_file <- function(fname,instrument='Cameca',suffix='.asc'){
    if (instrument=='Cameca' & suffix=='.asc'){
        out <- read_Cameca_asc(fname)
    } else if (instrument=='SHRIMP' & suffix=='.op'){
        out <- read_SHRIMP_op(fname)
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

subset_samples <- function(dat,prefix='Plesovice'){
    snames <- names(dat)
    matches <- grepl(prefix,snames)
    subset(dat,subset=matches)
}
