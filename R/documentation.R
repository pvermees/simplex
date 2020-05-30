#' @title example Cameca dataset
#' @description built-in Cameca U-Pb dataset containing Plesovice and
#'     Qinghu zircon
#' @name Cameca
#' @docType data
#' @examples
#' data(Cameca,package="simplex")
#' plot(Cameca,i=1)
NULL

#' @title example oxygen isotope dataset
#' @description built-in Cameca oxigen dataset containing NBS28,
#'     Qinghu zircon and ZBL10 data
#' @name oxygen
#' @docType data
#' @examples
#' data(oxygen,package="simplex")
#' plot(oxygen,i=1)
NULL

#' @title example SHRIMP dataset
#' @description built-in SHRIMP U-Pb dataset containing Temora and
#'     91500 zircon
#' @name SHRIMP
#' @docType data
#' @examples
#' data(SHRIMP,package="simplex")
#' plot(SHRIMP,i=1)
NULL

#' simplex: a SIMS data reduction package
#'
#' \code{simplex} is a data reduction package for Secondary Ion Mass
#' Spectrometry (SIMS), with functionality for both stable isotope
#' geochemistry and U-(Th)-Pb geochronology. The package accommodates
#' data from both Cameca and SHRIMP instruments. All mathematical
#' operations take place in logratio space (aka the simplex) to ensure
#' positivity and keep track of correlated uncertainties.
#'
#' A list of documented functions may be viewed by typing
#' \code{help(package='simplex')}.  Source code and further
#' instructions are provided at
#' \url{https://github.com/pvermees/simplex/}.
#' 
#' @name simplex
#' @docType package
#' @aliases simplex-package
"_PACKAGE"
#> [1] "_PACKAGE"
