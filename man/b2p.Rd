\name{b2p}
\alias{b2p}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
b2p(b, ag, tt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{ag}{
%%     ~~Describe \code{ag} here~~
}
  \item{tt}{
%%     ~~Describe \code{tt} here~~
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
function (b, ag, tt) 
{
    simplex <- c("204Pb", "206Pb", "207Pb", "238U", "238U16O2")
    nr <- nrow(tt)
    bm <- matrix(rep(b, nr), nrow = nr, byrow = TRUE)
    colnames(bm) <- simplex
    ebgt <- exp(bm + ag["g"] * tt)
    out <- ebgt/sum(ebgt)
    colnames(out) <- simplex
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
