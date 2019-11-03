\name{plot_timeresolved}
\alias{plot_timeresolved}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_timeresolved(samp, fit = FALSE, c64 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samp}{
%%     ~~Describe \code{samp} here~~
}
  \item{fit}{
%%     ~~Describe \code{fit} here~~
}
  \item{c64}{
%%     ~~Describe \code{c64} here~~
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
function (samp, fit = FALSE, c64 = NULL) 
{
    ions <- names(samp$dwelltime)
    np <- length(ions)
    nr <- ceiling(sqrt(np))
    nc <- ceiling(np/nr)
    par(mfrow = c(nr, nc), mai = c(0.4, 0.4, 0.1, 0.1))
    simplex <- c("204Pb", "206Pb", "207Pb", "238U", "238U16O2")
    if (fit) {
        X <- samp$time[, simplex]
        Y <- predict_cps(samp, c64 = c64)[, simplex]
    }
    for (ion in ions) {
        plot(samp$time[, ion], samp$cps[, ion], type = "p", xlab = "", 
            ylab = "")
        if (fit & ion \%in\% simplex) {
            lines(X[, ion], Y[, ion])
        }
        mtext(side = 1, text = "t", line = 2)
        mtext(side = 2, text = ion, line = 2)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
