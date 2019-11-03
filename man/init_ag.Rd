\name{init_ag}
\alias{init_ag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
init_ag(samp, c64 = NULL)
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
    t4 <- hours(samp$time[, "204Pb"])
    t6 <- hours(samp$time[, "206Pb"])
    c4 <- samp$counts[, "204Pb"]
    c6 <- samp$counts[, "206Pb"]
    fit6 <- glm(c6 ~ t6, family = poisson(link = "log"))
    b6 <- fit6$coef[1]
    g <- fit6$coef[2]
    b4 <- glm(c4 ~ 1 + offset(g * t4), family = poisson(link = "log"))$coef
    if (is.null(c64)) {
        a4 <- log(exp(b6 - b4) - 1)
    }
    else {
        d4 <- samp$dwelltime["204Pb"]
        d6 <- samp$dwelltime["206Pb"]
        a4 <- log(exp(b6 - b4) * (d4/d6)/c64 - 1)
    }
    out <- c(a4, g)
    names(out) <- c("a4", "g")
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
