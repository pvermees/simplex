\name{flatXYtable}
\alias{flatXYtable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
flatXYtable(XY)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XY}{
%%     ~~Describe \code{XY} here~~
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
function (XY) 
{
    snames <- names(XY)
    ns <- length(snames)
    xy <- matrix(0, ns, 5)
    colnames(xy) <- c("X", "sX", "Y", "sY", "rXY")
    for (i in 1:ns) {
        xy[i, "X"] <- XY[[i]]$x["X"]
        xy[i, "Y"] <- XY[[i]]$x["Y"]
        xy[i, "sX"] <- sqrt(XY[[i]]$cov["X", "X"])
        xy[i, "sY"] <- sqrt(XY[[i]]$cov["Y", "Y"])
        xy[i, "rXY"] <- XY[[i]]$cov["X", "Y"]/(xy[i, "sX"] * 
            xy[i, "sY"])
    }
    xy
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
