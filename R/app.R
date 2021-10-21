source("method.R")

presets <- function(method){
    defaultmethod(method)
}

freeformServer <- function(port=NULL) {
  appDir <- R.utils::getAbsolutePath("../inst/www")
  shinylight::slServer(host='0.0.0.0', port=port, appDir=appDir, daemonize=TRUE,
    interface=list(
      presets=presets
    )
  )
}

freeformServer(8000)
