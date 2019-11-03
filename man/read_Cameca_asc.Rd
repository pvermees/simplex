\name{read_Cameca_asc}
\alias{read_Cameca_asc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read_Cameca_asc(fname)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fname}{
%%     ~~Describe \code{fname} here~~
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
function (fname) 
{
    f <- file(fname)
    open(f)
    out <- list()
    while (length(line <- readLines(f, n = 1, warn = FALSE)) > 
        0) {
        if (grepl("ACQUISITION PARAMETERS", line)) {
            block <- readLines(f, n = 12, warn = FALSE)
            ions <- unlist(strsplit(block[2], split = "\t"))[-1]
            out$ions <- gsub(" ", "", ions)
            out$dwelltime <- as.numeric(unlist(strsplit(block[7], 
                split = "\t"))[-1])
            detector <- unlist(strsplit(block[12], split = "\t"))[-1]
            out$detector <- gsub(" ", "", detector)
            names(out$dwelltime) <- out$ions
            names(out$detector) <- out$ions
        }
        if (grepl("DETECTOR PARAMETERS", line)) {
            block <- readLines(f, n = 3, warn = FALSE)
            detectors <- NULL
            out$yield <- NULL
            out$background <- NULL
            while (TRUE) {
                line <- readLines(f, n = 1, warn = FALSE)
                if (grepl("CORRECTION FACTORS", line)) {
                  break
                }
                else if (nchar(line) > 0) {
                  detectorpars <- gsub(" ", "", unlist(strsplit(line, 
                    split = "\t")))
                  detectors <- c(detectors, detectorpars[1])
                  out$yield <- c(out$yield, as.numeric(detectorpars[2]))
                  out$background <- c(out$background, as.numeric(detectorpars[3]))
                }
            }
            names(out$yield) <- detectors
            names(out$background) <- detectors
        }
        if (grepl("RAW DATA", line)) {
            out$cps <- NULL
            junk <- readLines(f, n = 5, warn = FALSE)
            while ((line <- readLines(f, n = 1, warn = FALSE)) != 
                "") {
                dat <- as.numeric(unlist(strsplit(line, split = "\t")))[-c(1, 
                  2)]
                out$cps <- rbind(out$cps, dat)
            }
            colnames(out$cps) <- out$ions
            out$counts <- round(sweep(out$cps, MARGIN = 2, FUN = "*", 
                out$dwelltime))
        }
        if (grepl("PRIMARY INTENSITY", line)) {
            out$sbm <- NULL
            junk <- readLines(f, n = 5, warn = FALSE)
            while ((line <- readLines(f, n = 1, warn = FALSE)) != 
                "") {
                dat <- as.numeric(unlist(strsplit(line, split = "\t")))[-c(1, 
                  2)]
                out$sbm <- rbind(out$sbm, dat)
            }
            colnames(out$sbm) <- out$ions
        }
        if (grepl("TIMING", line)) {
            out$time <- NULL
            junk <- readLines(f, n = 5, warn = FALSE)
            while (length(line <- readLines(f, n = 1, warn = FALSE)) > 
                0) {
                dat <- as.numeric(unlist(strsplit(line, split = "\t")))[-c(1, 
                  2)]
                out$time <- rbind(out$time, dat)
            }
            colnames(out$time) <- out$ions
        }
    }
    close(f)
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
