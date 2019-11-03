\name{calibration}
\alias{calibration}
\title{
  Build a U-Pb calibration line
}

\description{ Regresses a line through
X=\eqn{^{238}}U\eqn{^{16}}O\eqn{_x}/\eqn{^{238}}U and
Y=(\eqn{^{206}}Pb-\eqn{^{204}}Pb[\eqn{^{206}}Pb/\eqn{^{204}}Pb]\eqn{_{\circ}})/\eqn{^{238}}U
}

\usage{ calibration(lr, dat, plot = TRUE, disp = TRUE, omit = NULL, ...)
}

\arguments{

  \item{lr}{
    a list of logratios and their covariance matrix, produced by \code{logratios}
  }
  
  \item{dat}{
    a data of class \code{simplex}
  }
  
  \item{plot}{
    set to \code{TRUE} if you want to plot the calibration diagram
  }
  
  \item{disp}{
    set to \code{TRUE} if you want to propagate the statistical
    uncertainty associated with the overdispersion of the standard data
  }
  
  \item{omit}{
    (vector of) integers marking the spots to omit from the regression analysis
  }
  
  \item{\dots}{
    optional arguments to be passed on to the generic \code{plot} function
  }
}

\value{
  A list with a vector \code{AB} containing the intercept and slope of
  the regression line; its covariance matrix; the MSWD of the fit and
  the degrees of freedom.
}

\seealso{
  \code{\link{calibrate},\link{logratios}}
}

\examples{
data(package="simplex")
stand <- subset_samples(dat=Cameca,prefix='Plesovice')
lr_stand <- logratios(dat=stand,c64=18.7)
fit <- calibration(lr_stand,dat=stand,plot=TRUE,disp=TRUE,omit=2)
}
