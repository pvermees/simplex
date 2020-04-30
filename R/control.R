.simplex <- new.env(parent = emptyenv())

.simplex$methods <- list()

set_method <- function(method,instrument,ions){
    .simplex$methods[[method]] <- list(instrument=instrument,ions=ions)
}

set_method(method='IGG-zircon',
           instrument='Cameca',
           ions=c('Zr90','Zr92','200.5','Zr94',
                  'Pb204','Pb206','Pb207','Pb208',
                  'U238','ThO2','UO2','270.1'))

set_method(method='IGG-monazite',
           instrument='Cameca',
           ions=c('La139','202.5','Pb204','Pb206',
                  'Pb207','Pb208','Th232','U238',
                  'ThO2','UO2'))

set_method(method='IGG-oxygen',
           instrument='Cameca',
           ions=c('O16','O17','O18'))

set_method(method='IGG-sulfur',
           instrument='Cameca',
           ions=c('S32','S33','33.96','S34','S36'))

set_method(method='GA-zircon',
           instrument='SHRIMP',
           ions=c('Zr2O','Pb204','bkg','Pb206','Pb207',
                  'Pb208','U238','ThO','UO','UO2'))

get_ions <- function(method){
    get_method(method=method,item='ions')
}
get_instrument <- function(method){
    get_method(method=method,item='instrument')
}
get_method <- function(method,item='ions'){
    if (method %in% names(.simplex$methods)) {
        return(.simplex$methods[[method]][[item]])
    } else {
        stop('Incorrect method.')
    }    
}
