\name{calplot}
\alias{calplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calplot(XY, fit, alpha = 0.05, sigdig = 2, omit = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XY}{
%%     ~~Describe \code{XY} here~~
}
  \item{fit}{
%%     ~~Describe \code{fit} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{sigdig}{
%%     ~~Describe \code{sigdig} here~~
}
  \item{omit}{
%%     ~~Describe \code{omit} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (XY, fit, alpha = 0.05, sigdig = 2, omit = NULL, ...) 
{
    xy <- flatXYtable(XY)
    isofit <- list(model = 1)
    isofit$fact <- stats::qt(1 - alpha/2, fit$df)
    isofit$n <- nrow(xy)
    isofit$mswd <- fit$mswd
    isofit$p.value <- as.numeric(1 - stats::pchisq(fit$mswd/fit$df, 
        fit$df))
    A <- fit$AB["A"]
    B <- fit$AB["B"]
    sA <- sqrt(fit$cov["A", "A"])
    sB <- sqrt(fit$cov["B", "B"])
    isofit$a <- c(A, sA, isofit$fact * sA, isofit$mswd * isofit$fact * 
        sA)
    isofit$b <- c(B, sB, isofit$fact * sB, isofit$mswd * isofit$fact * 
        sB)
    isofit$cov.ab <- fit$cov["A", "B"]
    xlab <- quote(""^238 * "U"^16 * "O"[2] * "/"^238 * "U")
    ylab <- quote(""^206 * "Pb*/"^238 * "U")
    scatterplot(xy, xlab = xlab, ylab = ylab, fit = isofit, omit = omit, 
        ...)
    graphics::title(isochrontitle(isofit, sigdig = sigdig), xlab = xlab, 
        ylab = ylab)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
