\name{raw_count_ratios}
\alias{raw_count_ratios}
\title{
  Time-resolved logratios
}
\description{
  Calculate the time-resolved logratios for a single spot.
}
\usage{
raw_count_ratios(samp)
}
\arguments{
  \item{samp}{
    One item in an object of class \code{simplex}.
  }
}
\value{
  A list with a vector of logratios and its covariance matrix.
}
\seealso{\code{\link{plot_timeresolved}}}
\examples{
data(package="simplex")
Lm <- raw_count_ratios(samp=Cameca[[1]])
}
