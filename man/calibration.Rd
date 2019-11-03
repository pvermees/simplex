\name{calibration}
\alias{calibration}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calibration(lr, dat, plot = TRUE, disp = TRUE, omit = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lr}{
%%     ~~Describe \code{lr} here~~
}
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
}
  \item{disp}{
%%     ~~Describe \code{disp} here~~
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
function (lr, dat, plot = TRUE, disp = TRUE, omit = NULL, ...) 
{
    XY <- lr2XY(lr = lr, dat = dat)
    fit <- yorkfit(XY, omit = omit)
    out <- list()
    out$AB <- fit$par
    out$df <- length(dat) - 2
    out$mswd <- 2 * fit$value/out$df
    if (disp) 
        d <- out$mswd
    else d <- 1
    out$cov <- d * solve(fit$hessian)
    labels <- c("A", "B")
    names(out$AB) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    if (plot) 
        calplot(XY, fit = out, omit = omit)
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
