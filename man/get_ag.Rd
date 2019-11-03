\name{get_ag}
\alias{get_ag}
\title{
  Prepare the bias correction
}
\description{
  Gets parameters \eqn{\alpha_4} and \eqn{\gamma} for the lag time bias
  correction of single collector SIMS data
}
\usage{
get_ag(samp, c64 = NULL)
}
\arguments{
  \item{samp}{
    One item in an object of class \code{simplex}.
  }
  \item{c64}{
    Composition (\eqn{^{206}}Pb/\eqn{^{204}}Pb ratio) of the common Pb.
  }
}
\value{
  A list containing a vector with parameters \eqn{\alpha_4} and
  \eqn{gamma}, and the corresponding covariance matrix.
}
\examples{
data(package="simplex")
samp <- Cameca[[1]]
Lm <- raw_count_ratios(samp=samp)
aLm <- avg_Lm(Lm)
ag <- get_ag(samp=samp)
out <- bias_correction(samp=samp,aLm=aLm,ag=ag)
}
