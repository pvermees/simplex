\name{simplex2isoplotr}
\alias{simplex2isoplotr}

\title{Save to \code{Isoplot}/\code{IsoplotR}}

\description{ Save the results in a flat table for further processing
  with \code{Isoplot} or \code{IsoplotR}.  }

\usage{ simplex2isoplotr(alr, format = 5) }

\arguments{
  
  \item{alr}{
    atomic logratio estimates (produced by \code{calibrate})
  }

  \item{format}{ \code{IsoplotR} input format. Currently only format 5
    is accommodated (\eqn{^{238}}U/\eqn{^{206}}Pb -
    \eqn{^{207}}Pb/\eqn{^{206}}Pb - \eqn{^{204}}Pb/\eqn{^{206}}Pb) } }

\value{
  a 9-column table of isotopic ratios, uncertainties and error correlations
}

\seealso{
  \code{\link{calibration},\link{calibrate}}
}

\examples{
data(package="simplex")
stand <- subset_samples(dat=Cameca,prefix='Plesovice')
lr_stand <- logratios(dat=stand,c64=18.7)
fit <- calibration(lr_stand,dat=stand,plot=FALSE,disp=TRUE,omit=2)
samp <- subset_samples(dat=Cameca,prefix='Qinghu')
lr_samp <- logratios(dat=samp)
logPbU <- calibrate(lr_samp,fit=fit,dat=samp,tst=c(337.13,0.18))
tab <- simplex2isoplotr(logPbU)
}
