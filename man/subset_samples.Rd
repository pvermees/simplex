\name{subset_samples}
\alias{subset_samples}

\title{
  Subset a dataset of class \code{simplex}
}

\description{
  Select a subset of samples aor standards from a \code{simplex} dataset.
}

\usage{
subset_samples(dat, prefix = "Plesovice")
}

\arguments{
  \item{dat}{
    A dataset to subset
  }
  
  \item{prefix}{
    Text string to match
  }
}

\value{
  an object with the same class as \code{dat}
}

\seealso{
  \code{\link{read_directory}}
}

\examples{
data(Cameca,package="simplex")
stand <- subset_samples(Cameca,prefix='Plesovice')
}