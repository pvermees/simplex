source("io.R")
source("method.R")

presets <- function(method){
    defaultmethod(method)
}

upload <- function(f,m){
    read_SHRIMP_op(f=textConnection(f),m=method(m))
}

freeformServer <- function(port=NULL) {
    appDir <- R.utils::getAbsolutePath("../inst/www")
    shinylight::slServer(host='0.0.0.0', port=port, appDir=appDir, daemonize=TRUE,
        interface=list(
          presets=presets,
          upload=upload
        )
    )
}

freeformServer(8000)
