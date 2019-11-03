\name{reshuffle}
\alias{reshuffle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
reshuffle(samps, c64 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samps}{
%%     ~~Describe \code{samps} here~~
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
function (samps, c64 = NULL) 
{
    out <- list(counts = list(), cps = list(), time = list(), 
        sbm = list(), dwelltime = NULL)
    snames <- names(samps)
    ns <- length(snames)
    ions <- names(samps[[1]]$dwelltime)
    for (i in 1:ns) {
        samp <- samps[[i]]
        for (ion in ions) {
            out$counts[[ion]] <- cbind(out$counts[[ion]], samp$counts[, 
                ion])
            out$cps[[ion]] <- cbind(out$cps[[ion]], samp$cps[, 
                ion] + 0.5/samp$dwelltime[ion])
            out$time[[ion]] <- cbind(out$time[[ion]], samp$time[, 
                ion])
            out$sbm[[ion]] <- cbind(out$sbm[[ion]], samp$sbm[, 
                ion])
        }
        out$dwelltime <- rbind(out$dwelltime, samp$dwelltime)
    }
    for (ion in ions) {
        colnames(out$counts[[ion]]) <- snames
        colnames(out$cps[[ion]]) <- snames
        colnames(out$time[[ion]]) <- snames
    }
    out$X <- log(out$cps[["238U16O2"]]) - log(out$cps[["238U"]])
    if (is.null(c64)) {
        out$Y <- log(out$cps[["206Pb"]]) - log(out$cps[["238U"]])
    }
    else {
        out$Y <- log(out$cps[["206Pb"]] - out$cps[["204Pb"]] * 
            c64) - log(out$cps[["238U"]])
    }
    colnames(out$X) <- snames
    colnames(out$Y) <- snames
    rownames(out$dwelltime) <- snames
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
