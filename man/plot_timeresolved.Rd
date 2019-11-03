\name{plot_timeresolved}
\alias{plot_timeresolved}
\title{
  Plot and fit raw time-resolved SIMS data.
}
\description{
  Shows the raw data for a single spot in a SIMS dataset.
}
\usage{
plot_timeresolved(samp, fit = FALSE, c64 = NULL)
}
\arguments{
  \item{samp}{
    One item of a \code{simplex} dataset.
  }
  \item{fit}{
    Set to \code{TRUE} to show the maximum likelihood fit on the plot.
  }
  \item{c64}{
    (optional) common Pb composition (\eqn{^206}Pb/\eqn{^204}Pb ratio).
  }
}
\value{
  Produces a multi-panel scatter plot.
}
\examples{
data(Cameca,package="simplex")
plot_timeresolved(Cameca[[1]],fit=TRUE)
}