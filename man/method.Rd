% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method.R
\name{method}
\alias{method}
\title{read or define SIMS acquisition details}
\usage{
method(m = "IGG-UPb", instrument, ions, num, den, bkg, description)
}
\arguments{
\item{m}{the name of a data acquisition protocol. Pre-defined
values include \code{'IGG-UPb'}, \code{'GA-UPb'},
\code{'IGG-UThPb'}, \code{'IGG-O'}, and \code{'IGG-S'}).}

\item{instrument}{one of either \code{SHRIMP} or \code{Cameca}.}

\item{ions}{vector of labels to be attached to the different ionic
masses that are visited during each sweep.}

\item{num}{numerator ions of the logratios to be processed in
subsequent data reduction steps.}

\item{den}{denominator ions of the logratios to be processed in
subsequent data reduction steps.}

\item{bkg}{either (1) a character string with the name of the
channel (e.g. \code{'bkg'}) of the background signal, (2)
\code{'nominal'} a flag indicating that nominal detector
backgrounds recorded in the input file are to be used (only
relevant to Cameca data), or (3) a numerical value representing
a nominal number of background counts (for SHRIMP) or counts
per second (for Cameca).}

\item{description}{text string with a description of the contents}
}
\value{
an object of class \code{method}
}
\description{
get or set a list of data acquisition properties
}
\examples{
fname <- system.file('SHRIMP.op',package='simplex')
m <- method(m="GA-UPb",
            instrument="SHRIMP",
            ions=c('Zr2O','Pb204','bkg','Pb206','Pb207',
                   'Pb208','U238','ThO','UO','UO2'),
            num=c('Pb204','Pb207','Pb206','UO'),
            den=c('Pb206','Pb206','U238','U238'),
            bkg='bkg',
            description="U-Pb dating at Geoscience Australia.")
shrimpdat <- read_data(fname,m=m)
plot(shrimpdat,i=1)
}
