d18O2ref <- function(d18O,err){
    out <- list()
    out$ref <- c(d18O/2,d18O)
    J <- matrix(0,2,1)
    J[1,1] <- 1/2
    J[2,1] <- 1
    out$cov <- J %*% (err^2) %*% t(J)
    out
}

standards_stable <- function(dat,prefix,invert=FALSE,ref,ref.cov){
    out <- list()
    out$ref <- ref
    out$ref.cov <- ref.cov
    out$x <- subset_samples(dat=dat,prefix=prefix,invert=invert)
    class(out) <- 'standard'
    out
}

calibration_stable <- function(stand,num=c('O17','O18'),den='O16',outliers=TRUE){
    ns <- length(stand$x)
    bclr <- list()
    for (i in 1:ns){
        spot <- mark_inliers(spot=stand$x[[i]],num=num,den=den)
        bclr[[i]] <- blankcor_multicol_spot(spot,num=num,den=den)
    }
    bclr
}
