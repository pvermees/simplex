\name{calibrate}
\alias{calibrate}

\title{
  Convert signal ratios to atomic ratios
}
\description{
  Compare the U-Pb ratios of a sample with that of the standard.
}
\usage{
calibrate(lr, fit, dat, PbUstand = NULL, tst = c(337.13, 0.18))
}

\arguments{

  \item{lr}{ a list of bias corrected logratios and their covariance
    matrices, produced by \code{logratios}.  }
  
  \item{fit}{ the output of \code{calibration}, containing the
    calibration curve of the standards data }

  \item{dat}{
    an object of class \code{simplex}
  }
  
  \item{PbUstand}{
    (optional) true \eqn{^{206}}Pb/\eqn{^{238}}U-ratio of the age standard
  }

  \item{tst}{
    list with the age and standard error of the age standard
  }
  
}

\value{ A list with the atomic \eqn{^{206}}Pb/\eqn{^{238}}U-ratios and
covariance matrix of all the samples.  }

\seealso{
  \code{\link{calibration},\link{logratios}}
}

\examples{
data(package="simplex")
stand <- subset_samples(dat=Cameca,prefix='Plesovice')
lr_stand <- logratios(dat=stand,c64=18.7)
fit <- calibration(lr_stand,dat=stand,plot=TRUE,disp=TRUE,omit=2)
samp <- subset_samples(dat=Cameca,prefix='Qinghu')
lr_samp <- logratios(dat=samp)
PbU <- calibrate(lr_samp,fit=fit,dat=samp,tst=c(337.13,0.18))
}
