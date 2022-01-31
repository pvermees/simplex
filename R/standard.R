#' @title define an isotopic reference standards
#' @description specify the isotopic composition of a reference
#'     material for SIMS data calibration.
#' @param preset (optional) text string. One of either
#' 
#' U-Pb age standards: \code{'Plesovice-t'}, \code{Qinghu-t},
#' \code{Penglai-t}, \code{'44069-t'}, \code{'Temora-t'},
#'
#' O-isotope standards: \code{'Plesovice-O'}, \code{'NBS28-O'},
#' \code{'Qinghu-O'}, \code{'Penglai-O'},
#'
#' or S-isotope standards: \code{'Sonora-S'}.
#' @param tst two element with the standard's age and its standard
#'     error.
#' @param measured logical. Only relevant for geochronological data.
#'     If \code{TRUE}, stores the composition of the standard as a
#'     measured isotopic composition. If \code{FALSE}, stores the
#'     composition as a mixture of common Pb and radiogenic Pb.
#' @param del Stable isotope composition of a standard in conventional
#'     delta notation, specified as a two-item list containing a
#'     vector \code{val} of \eqn{\delta}-values and its covariance
#'     matrix \code{cov}.
#' @param ref Logratio composition of a reference standard such as
#'     Vienna Standard Mean Ocean Water (VSMOW) or troilite. A
#'     two-item list containing a vector \code{val} of logratio-values
#'     and its covariance matrix \code{cov}.
#' @return an object of class \code{standard}
#' @examples
#' data(Cameca_UPb,package="simplex")
#' st <- standard(preset="Plesovice-t")
#' cal <- calibration(lr=Cameca_UPb,stand=st)
#' plot(cal)
#' @export
standard <- function(preset,tst,measured,del,ref){
    if (!missing(preset)){
        files <- system.file(c('tstand.csv','Ostand.csv','Sstand.csv'),
                             package='simplex')
        tpresets <- data.matrix(utils::read.csv(files[1],header=TRUE,
                                row.names='standard',check.names=FALSE))
        Opresets <- data.matrix(utils::read.csv(files[2],header=TRUE,
                                row.names='standard',check.names=FALSE))
        Spresets <- data.matrix(utils::read.csv(files[3],header=TRUE,
                                row.names='standard',check.names=FALSE))
        if (preset %in% rownames(tpresets)){
            geochron <- TRUE
            out <- age2stand(tst=tpresets[preset,])
        } else if (preset %in% rownames(Opresets)){
            geochron <- FALSE
            del <- list()
            del$val <- Opresets[preset,c(1,3)]
            del$cov <- diag(Opresets[preset,c(2,4)]^2)
            rownames(del$cov) <- colnames(del$cov) <- names(del$val)
            out <- del2stand(del,ref=VSMOW())
        } else if (preset %in% rownames(Spresets)){ # temporary value
            geochron <- FALSE
            del <- list()
            del$val <- Spresets[preset,c(1,3,5)]
            del$cov <- diag(Spresets[preset,c(2,4,6)]^2)
            rownames(del$cov) <- colnames(del$cov) <- names(del$val)
            out <- del2stand(del,ref=troilite())
        } else {
            stop("Invalid input to standard(...).")
        }
        out$preset <- preset
    } else if (!missing(tst)){
        geochron <- TRUE
        out <- age2stand(tst)
    } else if (!missing(del) & !missing(ref)){
        geochron <- FALSE
        out <- del2stand(del,ref)
    } else {
        stop('Illegal input to standard().')
    }
    if (geochron){
        if (missing(measured)) measured <- FALSE
        out$measured <- measured
        if (out$measured) {
            out$val <- out$val[-c(1:2)]
            out$cov <- out$cov[-c(1:2),-c(1:2)]
            if (!missing(preset)){
                if (preset=='Plesovice'){
                    out$val['Pb206/U238'] <- log(0.053707)
                    out$cov['Pb206/U238','Pb206/U238'] <- 0.02
                    out$cov['Pb206/U238','Pb208/Th232'] <- 0
                    out$cov['Pb208/Th232','Pb206/U238'] <- 0
                }
            }
        }
    } else {
        out$measured <- TRUE
    }
    out
}

age2stand <- function(tst){
    num <- c('Pb204','Pb204','Pb206','Pb208')
    den <- c('Pb206','Pb208','U238','Th232')
    ratios <- paste0(num,'/',den)
    common <- IsoplotR:::stacey.kramers(tst[1])
    D <- IsoplotR::mclean(tt=tst[1])
    val <- c(-log(common[,c('i64','i84')]),log(D$Pb206U238),log(D$Pb208Th232))
    J <- rbind(0,0,D$dPb206U238dt/D$Pb206U238,D$dPb208Th232dt/D$Pb208Th232)
    covmat <- J%*%(tst[2]^2)%*%t(J)
    names(val) <- ratios
    rownames(covmat) <- ratios
    colnames(covmat) <- ratios
    list(tst=tst,val=val,cov=covmat,measured=FALSE)
}

del2stand <- function(del,ref){
    keep <- (names(ref$val) %in% names(del$val))
    if (!any(keep)) stop('Standard isotopes must match reference.')
    out <- list()
    out$del <- list(val=del$val,cov=del$cov)
    out$ref <- list(preset=ref$preset,val=ref$val,cov=ref$cov)
    out$val <- log(1 + del$val/1000) + ref$val[keep]
    J <- diag(sum(keep))/(1000 + del$val)
    rownames(J) <- names(del$val)
    out$cov <- J %*% data.matrix(del$cov) %*% t(J)
    out
}

#' @title Vienna Standard Mean Ocean Water
#' @description returns the oxygen isotope composition of the
#'     eponymous stable isotope reference material
#' @return a two-item list containing a vector \code{val} of
#'     logratio-values and its covariance matrix \code{cov}.
#' @examples
#' data(Cameca_oxygen,package="simplex")
#' del <- list()
#' del$val <- 5.4
#' del$cov <- matrix(0.04,1,1)
#' names(del$val) <- colnames(del$cov) <- rownames(del$cov) <- 'O18/O16'
#' st <- standard(del=del,ref=VSMOW())
#' cal <- calibration(lr=Cameca_oxygen,stand=st)
#' plot(cal)
#' @export
VSMOW <- function(){
    ref(val=c(0.3799e-3,2.00520e-3),
        relerr=c(1.6e-3,0.43e-3)/c(0.3799,2.00520),
        preset='VSMOW',
        nms=c('O17/O16','O18/O16')
        )
}

#' @title troilite composition
#' @description returns the sulphur isotope composition of meteoritic
#'     troilite
#' @return a two-item list containing a vector \code{val} of
#'     logratio-values and its covariance matrix \code{cov}.
#' @examples
#' data(Cameca_sulphur,package="simplex")
#' del <- list()
#' del$val <- c(0.83,1.61,3.25)
#' del$cov <- diag(c(0.03,0.08,0.03)^2)
#' names(del$val) <- rownames(del$cov) <-
#' colnames(del$cov) <- c('S33/S32','S34/S32','S36/S32')
#' st <- standard(del=del,ref=troilite())
#' cal <- calibration(lr=Cameca_sulphur,stand=st)
#' plot(cal)
#' @export
troilite <- function(){
    S2346 <- c(126.948,22.6436,6515)
    ref(val=S2346,relerr=c(0.047,0.0020,20)/S2346,
        preset='troilite',
        nms=c('S33/S32','S34/S32','S36/S32'))
}

ref <- function(val,relerr,preset,nms){
    out <- list(val=val,cov=diag(relerr^2),preset=preset)
    names(out$val) <- colnames(out$cov) <- rownames(out$cov) <- nms
    out
}

skeletonstand <- function(lr,measured=TRUE){
    num <- lr$method$num
    den <- lr$method$den
    geochron <- ( all(c("Pb206","U238") %in% c(num,den)) |
                  all(c("Pb208","Th232") %in% c(num,den)) )
    num_is_element <- (element(num) %in% elements())
    den_is_element <- (element(den) %in% elements())
    ratios <- paste0(num,'/',den)
    if (geochron & measured){
        hasPb204 <- (num %in% "Pb204") | (den %in% "Pb204")
        selection <- ratios[num_is_element & den_is_element & !hasPb204]
    } else {
        selection <- ratios[num_is_element & den_is_element]
    }
    n2r <- length(selection)
    out <- list()
    out$measured <- measured
    out$val <- rep(0,n2r)
    out$cov <- matrix(0,n2r,n2r)
    names(out$val) <- selection
    rownames(out$cov) <- selection
    colnames(out$cov) <- selection
    out
}
