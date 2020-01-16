pars <- function(spot,oxide='UO2'){
    out <- list()
    out$t4 <- hours(spot$time[,'Pb204'])
    out$t6 <- hours(spot$time[,'Pb206'])
    out$t7 <- hours(spot$time[,'Pb207'])
    out$tU <- hours(spot$time[,'U238'])
    out$tO <- hours(spot$time[,oxide])
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
    if (spot$instrument=='Cameca'){
        bkg <- spot$background[spot$detector]
        names(bkg) <- names(spot$detector)
        nt <- length(out$t4)
        out$bc4 <- rep(bkg['Pb204'],nt)
        out$bc6 <- rep(bkg['Pb206'],nt)
        out$bc7 <- rep(bkg['Pb207'],nt)
        out$bcU <- rep(bkg['U238'],nt)
        out$bcO <- rep(bkg[oxide],nt)
        out$bd4 <- rep(1,nt)
        out$bd6 <- rep(1,nt)
        out$bd7 <- rep(1,nt)
        out$bdU <- rep(1,nt)
        out$bdO <- rep(1,nt)
        out$bt4 <- out$t4
        out$bt6 <- out$t6
        out$bt7 <- out$t7
        out$btU <- out$tU
        out$btO <- out$tO
        out$bn4 <- out$bc4*out$bt4
        out$bn6 <- out$bc6*out$bt6
        out$bn7 <- out$bc7*out$bt7
        out$bnU <- out$bcU*out$btU
        out$bnO <- out$bcO*out$btO
    } else {
        bt <- hours(spot$time[,'bkg'])
        bn <- spot$counts[,'bkg']
        bd <- spot$edt[,'bkg']
        bc <- bn/bd
        out$bc4 <- bc
        out$bc6 <- bc
        out$bc7 <- bc
        out$bcU <- bc
        out$bcO <- bc
        out$bd4 <- bd
        out$bd6 <- bd
        out$bd7 <- bd
        out$bdU <- bd
        out$bdO <- bd
        out$bt4 <- bt
        out$bt6 <- bt
        out$bt7 <- bt
        out$btU <- bt
        out$btO <- bt
        out$bn4 <- bn
        out$bn6 <- bn
        out$bn7 <- bn
        out$bnU <- bn
        out$bnO <- bn
    }
    out
}

hours <- function(tt){
    tt/3600
}
