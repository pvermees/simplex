\name{calibrate}
\alias{calibrate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calibrate(lr, fit, dat, PbUstand = NULL, tst = c(337.13, 0.18))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lr}{
%%     ~~Describe \code{lr} here~~
}
  \item{fit}{
%%     ~~Describe \code{fit} here~~
}
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{PbUstand}{
%%     ~~Describe \code{PbUstand} here~~
}
  \item{tst}{
%%     ~~Describe \code{tst} here~~
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
function (lr, fit, dat, PbUstand = NULL, tst = c(337.13, 0.18)) 
{
    XY <- lr2XY(lr = lr, dat = dat)
    xy <- flatXYtable(XY)
    ns <- nrow(xy)
    if (is.null(PbUstand)) 
        Pb6U8stand <- IsoplotR:::age_to_Pb206U238_ratio(tt = tst[1], 
            st = tst[2])
    Ytstand <- log(Pb6U8stand[1])
    out <- list()
    out$x <- xy[, "Y"] + Ytstand - fit$AB["A"] - fit$AB["B"] * 
        xy[, "X"]
    sYtstand <- Pb6U8stand[2]/Pb6U8stand[1]
    E <- matrix(0, 2 * ns + 3, 2 * ns + 3)
    i1 <- 1:ns
    i2 <- (ns + 1):(2 * ns)
    E[i1, i1] <- diag(xy[, "sX"]^2)
    E[i2, i2] <- diag(xy[, "sY"]^2)
    E[i1, i2] <- diag(xy[, "rXY"] * xy[, "sX"] * xy[, "sY"])
    E[i2, i1] <- E[i1, i2]
    E[(2 * ns + 1), (2 * ns + 1)] <- sYtstand^2
    E[(2 * ns) + (2:3), (2 * ns) + (2:3)] <- fit$cov
    J <- matrix(0, nrow = ns, ncol = 2 * ns + 3)
    rownames(J) <- rownames(xy)
    J[i1, i1] <- diag(-fit$AB["B"], ns, ns)
    J[i1, i2] <- diag(1, ns, ns)
    J[, 2 * ns + 1] <- 1
    J[, 2 * ns + 2] <- -1
    J[, 2 * ns + 3] <- -xy[, "X"]
    out$cov <- J \%*\% E \%*\% t(J)
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
