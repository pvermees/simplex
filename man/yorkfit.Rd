\name{yorkfit}
\alias{yorkfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
yorkfit(XY, omit = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XY}{
%%     ~~Describe \code{XY} here~~
}
  \item{omit}{
%%     ~~Describe \code{omit} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (XY, omit = NULL) 
{
    if (!is.null(omit)) {
        keep <- !((1:length(XY)) \%in\% omit)
        XY <- subset(XY, subset = keep)
    }
    snames <- names(XY)
    ns <- length(snames)
    X <- rep(0, ns)
    Y <- rep(0, ns)
    for (i in 1:ns) {
        X[i] <- XY[[i]]$x["X"]
        Y[i] <- XY[[i]]$x["Y"]
    }
    init <- lm(Y ~ X)$coef
    optim(init, misfit_york, method = "BFGS", XY = XY, hessian = TRUE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
