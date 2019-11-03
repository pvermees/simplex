\name{avg_Lm}
\alias{avg_Lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
avg_Lm(Lm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Lm}{
%%     ~~Describe \code{Lm} here~~
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
function (Lm) 
{
    np <- length(Lm[[1]]$x)
    init <- rep(0, np)
    nt <- length(Lm)
    for (i in 1:nt) {
        init <- init + Lm[[i]]$x/nt
    }
    lower <- init - rep(1, np)
    upper <- init + rep(1, np)
    fit <- optim(init, misfit_avg_Lm, Lm = Lm, lower = lower, 
        upper = upper, method = "L-BFGS-B", hessian = TRUE)
    out <- list()
    out$x <- fit$par
    out$cov <- solve(fit$hessian)
    labels <- c("L6m", "L7m", "LUm")
    names(out$x) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
