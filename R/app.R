setwd('/home/pvermees/Documents/Programming/R/simplex/')
source("R/io.R")
source("R/method.R")

presets <- function(method){
    if (method=='IGG-UPb'){
        load('data/Cameca_UPb.rda')
        out <- Cameca_UPb
    } else if (method=='IGG-UThPb'){
        load('data/Cameca_UThPb.rda')
        out <- Cameca_UThPb
    } else if (method=='IGG-O'){
        load('data/Cameca_oxygen.rda')
        out <- Cameca_oxygen
    } else if (method=='IGG-S'){
        load('data/Cameca_sulphur.rda')
        out <- Cameca_sulphur
    } else if (method=='GA-UPb'){
        load('data/SHRIMP_UPb.rda')
        out <- SHRIMP_UPb
    } else {
        out <- list()
        out$samples <- NULL
        out$method <- defaultmethod(method)
    }
    out
}

# f = list of two lists with blocks of text and corresponding filenames
# m = method (currently a string, will be modified to accept lists)
upload <- function(f,m){
    ntcs <- length(f$tcs)
    tcs <- list()
    for (i in 1:ntcs){
        tcs[[f$fns[[i]]]] <- textConnection(f$tcs[[i]])
    }
    read_data(f=tcs,m=m)
}

freeformServer <- function(port=NULL) {
    appDir <- R.utils::getAbsolutePath("inst/www")
    shinylight::slServer(host='0.0.0.0', port=port, appDir=appDir, daemonize=TRUE,
        interface=list(
          presets=presets,
          upload=upload
        )
    )
}

freeformServer(8000)
