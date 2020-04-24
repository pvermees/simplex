alpha <- function(spot,den){
    detector <- spot$detector[den]
    dtype <- spot$type[den]
    b <- spot$background[detector]
    d <- spot$counts[,den]
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
    n <- spot$counts[,num]
    d <- spot$counts[,den]
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

blankcor_multicol_spot <- function(spot,num,den,outliers=TRUE){
    nt <- nrow(spot$time)
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

inliers_multicol <- function(spot,num,den){
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
    inliers
}
