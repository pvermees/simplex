\name{ages}
\alias{ages}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ages(logPbU)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{logPbU}{
%%     ~~Describe \code{logPbU} here~~
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
function (logPbU) 
{
    out <- list()
    Pb6U8 <- exp(logPbU$x)
    JR <- diag(Pb6U8)
    E <- JR \%*\% logPbU$cov \%*\% t(JR)
    l38 <- IsoplotR::settings("lambda", "U238")[1]
    out$x <- log(1 + Pb6U8)/l38
    Jt <- diag(1/(l38 * (1 + Pb6U8)))
    out$cov <- Jt \%*\% E \%*\% t(Jt)
    labels <- names(logPbU$x)
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
