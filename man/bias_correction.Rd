\name{bias_correction}
\alias{bias_correction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bias_correction(samp, aLm, ag, c64 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samp}{
%%     ~~Describe \code{samp} here~~
}
  \item{aLm}{
%%     ~~Describe \code{aLm} here~~
}
  \item{ag}{
%%     ~~Describe \code{ag} here~~
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
function (samp, aLm, ag, c64 = NULL) 
{
    simplex <- c("204Pb", "206Pb", "207Pb", "238U", "238U16O2")
    tt <- hours(samp$time[, simplex])
    dt6 <- mean(tt[, "206Pb"] - tt[, "238U"])
    dt7 <- mean(tt[, "207Pb"] - tt[, "206Pb"])
    dtU <- mean(tt[, "238U16O2"] - tt[, "238U"])
    L6c <- aLm$x["L6m"] - ag$x["g"] * dt6
    L7c <- aLm$x["L7m"] - ag$x["g"] * dt7
    LUc <- aLm$x["LUm"] - ag$x["g"] * dtU
    if (is.null(c64)) {
        L4c <- log(exp(ag$x["a4"]) + 1)
    }
    else {
        d4 <- samp$dwelltime["204Pb"]
        d6 <- samp$dwelltime["206Pb"]
        L4c <- log(exp(ag$x["a4"]) + 1) + log(c64 * d6/d4)
    }
    out <- list()
    out$x <- c(L4c, L6c, L7c, LUc)
    labels <- c("L4c", "L6c", "L7c", "LUc")
    names(out$x) <- labels
    E <- matrix(0, 5, 5)
    E[1:3, 1:3] <- aLm$cov
    E[4:5, 4:5] <- ag$cov
    J <- matrix(0, 4, 5)
    rownames(J) <- labels
    colnames(J) <- c("L6m", "L7m", "LUm", "a4", "g")
    J["L6c", "L6m"] <- 1
    J["L7c", "L7m"] <- 1
    J["LUc", "LUm"] <- 1
    J["L4c", "a4"] <- exp(ag$x["a4"])/(exp(ag$x["a4"]) + 1)
    J["L6c", "g"] <- -dt6
    J["L7c", "g"] <- -dt7
    J["LUc", "g"] <- -dtU
    out$cov <- J \%*\% E \%*\% t(J)
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
