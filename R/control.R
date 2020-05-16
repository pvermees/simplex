.simplex <- new.env(parent = emptyenv())

.simplex$methods <- list()

set_method <- function(method,type,instrument,ions,nominalblank){
    .simplex$methods[[method]] <- list(instrument=instrument,type=type,
                                       ions=ions,nominalblank=nominalblank)
}

set_method(method='IGG-zircon',
           type='U-Pb',
           instrument='Cameca',
           ions=c('Zr90','Zr92','200.5','Zr94',
                  'Pb204','Pb206','Pb207','Pb208',
                  'U238','ThO2','UO2','270.1'),
           nominalblank=TRUE)

set_method(method='IGG-monazite',
           type='U-Th-Pb',
           instrument='Cameca',
           ions=c('La139','202.5','Pb204','Pb206',
                  'Pb207','Pb208','Th232','U238',
                  'ThO2','UO2'),
           nominalblank=TRUE)

set_method(method='IGG-oxygen',
           type='oxygen',
           instrument='Cameca',
           ions=c('O16','O17','O18'),
           nominalblank=TRUE)

set_method(method='IGG-sulfur',
           type='sulfur',
           instrument='Cameca',
           ions=c('S32','S33','33.96','S34','S36'),
           nominalblank=TRUE)

set_method(method='GA-zircon',
           instrument='SHRIMP',
           type='U-Pb',
           ions=c('Zr2O','Pb204','bkg','Pb206','Pb207',
                  'Pb208','U238','ThO','UO','UO2'),
           nominalblank=FALSE)

get_ions <- function(method){
    get_method(method=method,item='ions')
}
get_instrument <- function(method){
    get_method(method=method,item='instrument')
}
get_type <- function(method){
    get_method(method=method,item='type')
}
nominalblank <- function(method){
    get_method(method=method,item='nominalblank')
}
get_method <- function(method,item='ions'){
    if (method %in% names(.simplex$methods)) {
        return(.simplex$methods[[method]][[item]])
    } else {
        stop('Incorrect method.')
    }    
}
