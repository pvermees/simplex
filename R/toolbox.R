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
    out$d4 <- spot$edt[,'Pb204']
    out$d6 <- spot$edt[,'Pb206']
    out$d7 <- spot$edt[,'Pb207']
    out$dU <- spot$edt[,'U238']
    out$dO <- spot$edt[,oxide]
    out$c4 <- out$n4/out$d4
    out$c6 <- out$n6/out$d6
    out$c7 <- out$n7/out$d7
    out$cU <- out$nU/out$dU
    out$cO <- out$nO/out$dO
    out
}
