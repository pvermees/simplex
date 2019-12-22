pars <- function(spot,oxide='UO2'){
    out <- list()
    out$t4 <- spot$time[,'Pb204']
    out$t6 <- spot$time[,'Pb206']
    out$t7 <- spot$time[,'Pb207']
    out$tU <- spot$time[,'U238']
    out$tO <- spot$time[,oxide]
    out$n4 <- spot$counts[,'Pb204']
    out$n6 <- spot$counts[,'Pb206']
    out$n7 <- spot$counts[,'Pb207']
    out$nU <- spot$counts[,'U238']
    out$nO <- spot$counts[,oxide]
    out$d4 <- spot$dwelltime['Pb204']
    out$d6 <- spot$dwelltime['Pb206']
    out$d7 <- spot$dwelltime['Pb207']
    out$dU <- spot$dwelltime['U238']
    out$dO <- spot$dwelltime[oxide]
    # dead time correction:
    dt <- spot$deadtime[spot$detector]
    names(dt) <- names(spot$detector)
    out$c4 <- out$n4/(out$d4-dt['Pb204']*out$n4*1e-9)
    out$c6 <- out$n6/(out$d6-dt['Pb206']*out$n6*1e-9)
    out$c7 <- out$n7/(out$d7-dt['Pb207']*out$n7*1e-9)
    out$cU <- out$nU/(out$dU-dt['U238']*out$nU*1e-9)
    out$cO <- out$nO/(out$dO-dt[oxide]*out$nO*1e-9)
    out
}
