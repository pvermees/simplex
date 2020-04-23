get_ions <- function(format){
    if (format == 'IGG-zircon'){
        out <- c('Zr90','Zr92','200.5','Zr94',
                 'Pb204','Pb206','Pb207','Pb208',
                 'U238','ThO2','UO2','270.1')
    } else if (format == 'IGG-monazite'){
        out <- c('La139','202.5','Pb204','Pb206',
                 'Pb207','Pb208','Th232','U238',
                 'ThO2','UO2')
    } else if (format == 'IGG-oxygen'){
        out <- c('O16','O17','O18')
    } else if (format == 'IGG-sulfur'){
        out <- c('S32','S33','33.96','S34','S36')
    } else if (format == 'GA-zircon'){
        out <- c('Zr2O','Pb204','bkg','Pb206','Pb207',
                 'Pb208','U238','ThO','UO','UO2')
    } else {
        stop('Incorrect format.')
    }
    out
}
