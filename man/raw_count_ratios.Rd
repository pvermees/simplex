\name{raw_count_ratios}
\alias{raw_count_ratios}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
raw_count_ratios(samp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samp}{
%%     ~~Describe \code{samp} here~~
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
function (samp) 
{
    counts <- samp$counts[, c("206Pb", "207Pb", "238U", "238U16O2")]
    lc <- log(counts)
    nt <- nrow(lc)
    n6 <- counts[, "206Pb"]
    n7 <- counts[, "207Pb"]
    nU <- counts[, "238U"]
    nUO <- counts[, "238U16O2"]
    L6m <- lc[, "206Pb"] - lc[, "238U"]
    L7m <- lc[, "207Pb"] - lc[, "206Pb"]
    LUm <- lc[, "238U16O2"] - lc[, "238U"]
    rs <- rowSums(counts)
    D <- exp(L6m) + exp(L7m + L6m) + exp(LUm) + 1
    dDdL6m <- exp(L6m) + exp(L7m + L6m)
    dDdL7m <- exp(L7m + L6m)
    dDdLUm <- exp(LUm)
    d2DdL6m2 <- exp(L6m) + exp(L7m + L6m)
    d2DdL7m2 <- exp(L7m + L6m)
    d2DdLUm2 <- exp(LUm)
    d2DdL6mdL7m <- exp(L7m + L6m)
    d2DdL7mdL6m <- exp(L7m + L6m)
    d2DdL6mdLUm <- 0
    d2DdL7mdLUm <- 0
    d2DdLUmdL6m <- 0
    d2DdLUmdL7m <- 0
    d2LLdL6m2 <- (rs/D^2) * (dDdL6m)^2 - (rs/D) * d2DdL6m2
    d2LLdL6mdL7m <- (rs/D^2) * dDdL6m * dDdL7m - (rs/D) * d2DdL6mdL7m
    d2LLdL6mdLUm <- (rs/D^2) * dDdL6m * dDdLUm - (rs/D) * d2DdL6mdLUm
    d2LLdL7m2 <- (rs/D^2) * (dDdL7m)^2 - (rs/D) * d2DdL7m2
    d2LLdL7mdL6m <- (rs/D^2) * dDdL7m * dDdL6m - (rs/D) * d2DdL6mdL7m
    d2LLdL7mdLUm <- (rs/D^2) * dDdL7m * dDdLUm - (rs/D) * d2DdL7mdLUm
    d2LLdLUm2 <- (rs/D^2) * (dDdLUm)^2 - (rs/D) * d2DdLUm2
    d2LLdLUmdL6m <- (rs/D^2) * dDdLUm * dDdL6m - (rs/D) * d2DdLUmdL6m
    d2LLdLUmdL7m <- (rs/D^2) * dDdLUm * dDdL7m - (rs/D) * d2DdLUmdL7m
    out <- list()
    labels <- c("L6m", "L7m", "LUm")
    for (i in 1:nt) {
        out[[i]] <- list()
        out[[i]]$x <- c(L6m[i], L7m[i], LUm[i])
        H <- matrix(0, 3, 3)
        H[1, 1] <- d2LLdL6m2[i]
        H[1, 2] <- d2LLdL6mdL7m[i]
        H[1, 3] <- d2LLdL6mdLUm[i]
        H[2, 1] <- d2LLdL7mdL6m[i]
        H[2, 2] <- d2LLdL7m2[i]
        H[2, 3] <- d2LLdL7mdLUm[i]
        H[3, 1] <- d2LLdLUmdL6m[i]
        H[3, 2] <- d2LLdLUmdL7m[i]
        H[3, 3] <- d2LLdLUm2[i]
        out[[i]]$H <- H
        names(out[[i]]$x) <- labels
        colnames(out[[i]]$H) <- labels
        rownames(out[[i]]$H) <- labels
    }
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
