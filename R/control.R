.simplex <- new.env(parent = emptyenv())

.simplex$ions <- list(
    'IGG-zircon' = c('Zr90','Zr92','200.5','Zr94',
                     'Pb204','Pb206','Pb207','Pb208',
                     'U238','ThO2','UO2','270.1'),
    'IGG-monazite' = c('La139','202.5','Pb204','Pb206',
                       'Pb207','Pb208','Th232','U238',
                       'ThO2','UO2'),
    'IGG-oxygen' =  c('O16','O17','O18'),
    'IGG-sulfur' = c('S32','S33','33.96','S34','S36'),
    'GA-zircon' = c('Zr2O','Pb204','bkg','Pb206','Pb207',
                    'Pb208','U238','ThO','UO','UO2')
)

set_method <- function(method,ions){
    .simplex$ions[[method]] <- ions
}

get_ions <- function(method){
    if (method %in% names(.simplex$ions)) {
        return(.simplex$ions[[method]])
    } else {
        stop('Incorrect method.')
    }
}
