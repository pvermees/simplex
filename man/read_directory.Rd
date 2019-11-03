\name{read_directory}
\alias{read_directory}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read_directory(dname, instrument = "Cameca", suffix = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dname}{
%%     ~~Describe \code{dname} here~~
}
  \item{instrument}{
%%     ~~Describe \code{instrument} here~~
}
  \item{suffix}{
%%     ~~Describe \code{suffix} here~~
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
function (dname, instrument = "Cameca", suffix = NULL) 
{
    if (instrument == "Cameca") {
        if (is.null(suffix)) 
            suffix <- ".asc"
    }
    else if (instrument == "SHRIMP") {
        if (is.null(suffix)) 
            suffix <- ".op"
    }
    else {
        stop("Unsupported instrument")
    }
    out <- list()
    class(out) <- instrument
    fnames <- list.files(dname, pattern = suffix)
    nf <- length(fnames)
    for (i in 1:nf) {
        fname <- fnames[i]
        sname <- tools::file_path_sans_ext(fname)
        out[[sname]] <- read_file(paste0(dname, fname), instrument = instrument, 
            suffix = suffix)
    }
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
