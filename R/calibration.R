calplot <- function(b){
    B <- flatten_beta(b)
    X <- B[,'UO2/U238']
    Y <- B[,'Pb206/U238']
    plot(X,Y)
}

flatten_beta <- function(b){
    snames <- names(b)
    ns <- length(snames)
    nc <- ncol(b[[1]])
    out <- matrix(0,ns,nc)
    colnames(out) <- colnames(b[[1]])
    rownames(out) <- snames
    for (sname in snames){
        out[sname,] <- b[[sname]]['b0',]
    }
    out
}
