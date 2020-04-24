alpha <- function(spot,den){
    detector <- spot$detector[den]
    dtype <- spot$type[den]
    b <- spot$background[detector]
    d <- spot$counts[inliers(spot),den]
    if (dtype=='Fc'){
        out <- log(1-b/d)
    } else if (dtype=='Em') {
        stop('Multicollection for electron multipliers not implemented yet.')
    } else {
        stop('Illegal detector type.')
    }
    out
}

beta <- function(spot,num,den){
    dnum <- spot$detector[num]
    dden <- spot$detector[den]
    dntype <- spot$type[num]
    ddtype <- spot$type[den]
    i <- inliers(spot)
    n <- spot$counts[i,num]
    d <- spot$counts[i,den]
    if (dntype=='Fc' & ddtype=='Fc'){
        out <- log(n) - log(d)
    } else if (dntype=='Em' & ddtype=='Fc') {
        stop('Mixed Fc and Em dectors not implemented yet.')
    } else if (dntype=='Fc' & ddtype=='Em') {
        stop('Mixed Fc and Em dectors not implemented yet.')
    } else if (dntype=='Em' & ddtype=='Em'){
        stop('Multicollection for electron multipliers not implemented yet.')
    } else {
        stop('Illegal detector type.')
    }
    out
}

blankcor_multicol_spot <- function(spot,num,den){
    nt <- length(inliers(spot))
    nn <- length(num)
    b <- matrix(0,nt,nn)
    colnames(b) <- num
    for (n in num){
        b[,n] <- beta(spot,n,den) + alpha(spot,n) - alpha(spot,den)
    }
    out <- list()
    out$b <- colMeans(b)
    out$cov <- cov(b)/nt
    out
}

inliers_multicol <- function(spot,num,den,plot=FALSE){
    lr <- sweep(log(spot$counts[,num,drop=FALSE]),1,log(spot$counts[,den]))
    nt <- nrow(lr)
    inliers <- 1:nt
    df <- length(num)
    while (TRUE){
        LR <- lr[inliers,,drop=FALSE]
        MD <- mahalanobis(x=LR,center=colMeans(LR),cov=cov(LR))
        imaxMD <- which.max(MD)
        cutoff <- qchisq(1-0.5/nt,df=df)
        if (MD[imaxMD] > cutoff){
            inliers <- inliers[-imaxMD]
        } else {
            break
        }
    }
    if (plot){
        pch <- rep(4,nt)
        pch[inliers] <- 1
        plot(lr,pch=pch)
    }
    inliers
}

mark_inliers <- function(spot,num,den){
    out <- spot
    out$inliers <- inliers_multicol(spot=spot,num=num,den=den)
    out
}

inliers <- function(spot){
    if ('inliers' %in% names(spot)){
        out <- spot$inliers
    } else {
        out <- 1:nrow(spot$time)
    }
    out
}
