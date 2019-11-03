\name{logratios}
\alias{logratios}

\title{
  Logratios for all the spots in a dataset
}

\description{
  Create a list of bias corrected logratios and their covariance matrices.
}

\usage{ logratios(dat, c64 = NULL) }

\arguments{
  
  \item{dat}{
    an object of class \code{simplex}
  }

  \item{c64}{
    composition (\eqn{^{206}}Pb/\eqn{^{204}}Pb-ratio) of the common Pb.
  }
}

\value{
  a list of bias corrected logratios and their covariance matrices.
}

\seealso{
  \code{\link{calibration},\link{calibrate}}
}

\examples{
data(package="simplex")
stand <- subset_samples(dat=Cameca,prefix='Plesovice')
lr_stand <- logratios(dat=stand,c64=18.7)
}
