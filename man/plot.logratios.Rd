% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logratios.R
\name{plot.logratios}
\alias{plot.logratios}
\title{plot logratios}
\usage{
\method{plot}{logratios}(x, sname = NULL, i = 1, ratios = FALSE, ...)
}
\arguments{
\item{x}{an object of class \code{logratios}}

\item{sname}{the sample name to be shown}

\item{i}{the sample number to be shown}

\item{ratios}{logical. If \code{FALSE}, plots the raw signals
versus time. If \code{TRUE}, plots the ratios against time.
Both plots show the fitted values as a solid line. Note that,
for single collector datasets, the numerator and denominator of
the measured ratios correspond to different times.}

\item{...}{optional arguments to be passed on to the generic
\code{plot} function.}
}
\description{
plot time resolved logratio data
}
\examples{
\dontrun{
data('SHRIMP_UPb',package='simplex')
dc <- drift(x=SHRIMP_UPb)
lr <- logratios(dc)
plot(lr,i=1,ratios=TRUE)
}
}
