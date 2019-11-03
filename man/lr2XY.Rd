\name{lr2XY}
\alias{lr2XY}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lr2XY(lr, dat, c64 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lr}{
%%     ~~Describe \code{lr} here~~
}
  \item{dat}{
%%     ~~Describe \code{dat} here~~
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
function (lr, dat, c64 = NULL) 
{
    snames <- names(dat)
    ns <- length(snames)
    out <- list()
    labels <- c("X", "Y")
    calions <- c("L4c", "L6c", "LUc")
    J <- matrix(0, 2, 3)
    rownames(J) <- labels
    colnames(J) <- calions
    J["X", "LUc"] <- 1
    J["Y", "L6c"] <- 1
    for (sname in snames) {
        d4 <- dat[[sname]]$dwelltime["204Pb"]
        d6 <- dat[[sname]]$dwelltime["206Pb"]
        dU <- dat[[sname]]$dwelltime["238U"]
        dUO <- dat[[sname]]$dwelltime["238U16O2"]
        L4c <- lr[[sname]]$x["L4c"]
        L6c <- lr[[sname]]$x["L6c"]
        LUc <- lr[[sname]]$x["LUc"]
        X <- LUc + log(dU/dUO)
        if (is.null(c64)) {
            Y <- L6c + log(dU/d6)
        }
        else {
            Y <- log(dU) + log(dU/d6) + log(1 - (c64 * d6)/(d4 * 
                exp(L4c)))
            J["Y", "L4c"] <- c64 * d6/(c64 * d6 + exp(L4c) * 
                d4)
        }
        out[[sname]] <- list()
        out[[sname]]$x <- c(X, Y)
        names(out[[sname]]$x) <- labels
        out[[sname]]$cov <- J \%*\% lr[[sname]]$cov[calions, calions] \%*\% 
            t(J)
        out[[sname]]$omega <- solve(out[[sname]]$cov)
    }
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
