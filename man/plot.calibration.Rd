% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibration.R
\name{plot.calibration}
\alias{plot.calibration}
\title{plot calibration data}
\usage{
\method{plot}{calibration}(x, option = 1, snames = NULL, i = NULL, ...)
}
\arguments{
\item{x}{an object of class \code{logratios}}

\item{option}{if \code{option=1}, plots the best fit line through
U-Pb and Th-Pb data. If \code{option=2}, adds the time-resolved
raw data to the plot. If \code{option=3}, marks the first and
last measurement by black and white circles, respectively.}

\item{snames}{the sample names to be shown}

\item{i}{the sample number to be shown}

\item{...}{optional arguments to be passed on to the generic
\code{plot} function.}
}
\description{
shows the calibration data on a logratio plot.
}
\examples{
\dontrun{
data('SHRIMP',package='simplex')
st <- standard(preset='Temora',prefix=TEM)
dc <- drift(x=SHRIMP)
lr <- logratios(x=dc)
cal <- calibration(lr=lr,stand=st)
plot(cal,option=3)
}
}
