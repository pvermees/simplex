\name{predict_cps}
\alias{predict_cps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict_cps(samp, c64 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
function (samp, c64 = NULL) 
{
    simplex <- c("204Pb", "206Pb", "207Pb", "238U", "238U16O2")
    Lm <- raw_count_ratios(samp = samp)
    aLm <- avg_Lm(Lm)
    ag <- get_ag(samp = samp, c64 = c64)
    if (is.null(c64)) {
        L4m <- log(exp(ag$x["a4"]) + 1)
    }
    else {
        d4 <- samp$dwelltime["204Pb"]
        d6 <- samp$dwelltime["206Pb"]
        L4m <- log(exp(ag$x["a4"]) + 1) + log(c64 * d6/d4)
    }
    L6m <- aLm$x["L6m"]
    L7m <- aLm$x["L7m"]
    LUm <- aLm$x["LUm"]
    den <- exp(L6m - L4m) + exp(L6m) + exp(L7m + L6m) + exp(LUm) + 
        1
    p4 <- exp(L6m - L4m)/den
    p6 <- exp(L6m)/den
    p7 <- exp(L7m + L6m)/den
    pU <- 1/den
    pUO <- exp(LUm)/den
    p <- c(p4, p6, p7, pU, pUO)
    rs <- rowSums(samp$counts[, simplex])
    dt <- samp$dwelltime[simplex]
    counts <- matrix(rs, ncol = 1) \%*\% matrix(p, nrow = 1)
    cps <- sweep(counts, MARGIN = 2, FUN = "/", dt)
    colnames(cps) <- simplex
    cps
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
