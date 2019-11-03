\name{misfit_york}
\alias{misfit_york}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
misfit_york(AB, XY = XY)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{AB}{
%%     ~~Describe \code{AB} here~~
}
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
function (AB, XY = XY) 
{
    A <- AB[1]
    B <- AB[2]
    SS <- 0
    snames <- names(XY)
    D <- matrix(0, 1, 2)
    for (sname in snames) {
        X <- XY[[sname]]$x["X"]
        Y <- XY[[sname]]$x["Y"]
        O <- XY[[sname]]$omega
        CC <- Y - A - B * X
        C1 <- O[1, 1] + O[1, 2] * B + O[2, 1] * B + O[2, 2] * 
            B^2
        C2 <- 2 * (O[1, 2] + O[2, 2] * B) * CC
        C3 <- O[2, 2] * CC^2
        K <- -C2/(2 * C1)
        SS <- SS + C1 * K^2 + C2 * K + C3
    }
    SS/2
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
