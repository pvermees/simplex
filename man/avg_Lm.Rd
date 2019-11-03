\name{avg_Lm}
\alias{avg_Lm}
\title{
  Average logratios
}
\description{
  Takes time resolved logratio estimates as input and produces a single
  spot logratio vector and covariance matrix as output.
}
\usage{
avg_Lm(Lm)
}

\arguments{
  \item{Lm}{
    A time resolved list of logratio vectors and Hessian matrices.
  }
}
\details{
  Uses a weighted mean algorithm with correlated uncertainties.
}
\value{
  A list with a single vector of logratios and its covariance matrix.
}
\examples{
data(Cameca,package='simplex')
Lm <- raw_count_ratios(samp=Cameca[[1]])
aLm <- avg_Lm(Lm)
}