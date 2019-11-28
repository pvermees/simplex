data2table <- function(x){
    UPb <- logratios2ratios(x)
    ns <- length(x$names)
    nr <- length(x$logratios)
    out <- matrix(0,ns,3*nr)
    rownames(out) <- x$names
    colnames(out) <- c('U238Pb206','s[U238Pb206]',
                       'Pb207Pb206','s[Pb207Pb206]',
                       'Pb204Pb206','s[Pb204Pb206]',
                       'rXY','rXZ','rYZ')
    cormat <- cov2cor(UPb$cov)
    for (i in 1:ns){
        j <- (i-1)*nr+(1:nr)
        out[i,c('U238Pb206','Pb207Pb206',
                'Pb204Pb206')] <- UPb$x[j]
        out[i,c('s[U238Pb206]','s[Pb207Pb206]',
                's[Pb204Pb206]')] <- sqrt(diag(UPb$cov))[j]
        out[i,'rXY'] <- cormat[j[1],j[2]]
        out[i,'rXZ'] <- cormat[j[1],j[3]]
        out[i,'rYZ'] <- cormat[j[2],j[3]]
    }
    out
}
