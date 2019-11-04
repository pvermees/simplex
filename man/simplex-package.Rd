\name{simplex-package}
\alias{simplex-package}
\alias{simplex}
\docType{package}
\title{
\packageTitle{simplex}
}
\description{
\packageDescription{simplex}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{simplex}
\packageIndices{simplex}
}
\author{
\packageAuthor{simplex}

Maintainer: \packageMaintainer{simplex}
}

\examples{
data(package="simplex")
stand <- subset_samples(dat=Cameca,prefix='Plesovice')
lr_stand <- logratios(dat=stand,c64=18.7)
fit <- calibration(lr_stand,dat=stand,plot=TRUE,disp=TRUE,omit=2)
samp <- subset_samples(dat=Cameca,prefix='Qinghu')
lr_samp <- logratios(dat=samp)
PbU <- calibrate(lr_samp,fit=fit,dat=samp,tst=c(337.13,0.18))
tab <- simplex2isoplotr(logPbU)
}