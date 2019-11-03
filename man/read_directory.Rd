\name{read_directory}
\alias{read_directory}
\title{
  Read data directory
}
\description{
  Read all the input files in a data directory
}
\usage{
read_directory(dname, instrument = "Cameca", suffix = NULL)
}
\arguments{
  
  \item{dname}{
    path to the input directory
  }
  
  \item{instrument}{
    text string with the type of ICP-MS. Currently only
    \code{'Cameca'}.
  }
  
  \item{suffix}{
    (optional) file extension of the input files.
  }
}
\value{
  An object of class \code{simplex}, i.e. a list of lists
  containing the following items: \code{ions}, \code{ions},
  \code{dwelltime}, \code{detector}, \code{yield}, \code{background},
  \code{cps}, \code{counts}, \code{sbm}, and \code{time}.
}
\examples{
datadir <- system.file(package="simplex")
dat <- read_directory(paste0(datadir,'/'),instrument='Cameca',suffix='.asc')
}