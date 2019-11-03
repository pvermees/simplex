\name{ag2b}
\alias{ag2b}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ag2b(ag, samp, c64 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ag}{
%%     ~~Describe \code{ag} here~~
}
  \item{samp}{
%%     ~~Describe \code{samp} here~~
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
function (ag, samp, c64 = NULL) 
{
    simplex <- c("204Pb", "206Pb", "207Pb", "238U", "238U16O2")
    tt <- hours(samp$time[, simplex])
    cc <- samp$counts[, simplex]
    dt <- samp$dwelltime[simplex]
    fn6 <- cc[, "206Pb"] ~ 1 + offset(ag["g"] * tt[, "206Pb"])
    fn7 <- cc[, "207Pb"] ~ 1 + offset(ag["g"] * tt[, "207Pb"])
    fnU <- cc[, "238U"] ~ 1 + offset(ag["g"] * tt[, "238U"])
    fnUO <- cc[, "238U16O2"] ~ 1 + offset(ag["g"] * tt[, "238U16O2"])
    b6 <- glm(fn6, family = poisson(link = "log"))$coef
    b7 <- glm(fn7, family = poisson(link = "log"))$coef
    bU <- glm(fnU, family = poisson(link = "log"))$coef
    bUO <- glm(fnUO, family = poisson(link = "log"))$coef
    if (is.null(c64)) {
        b4 <- b6 - log(exp(ag["a4"]) + 1)
    }
    else {
        b4 <- b6 - log(c64 * dt["206Pb"]/dt["204Pb"]) - log(exp(ag["a4"]) + 
            1)
    }
    out <- c(b4, b6, b7, bU, bUO)
    names(out) <- simplex
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
