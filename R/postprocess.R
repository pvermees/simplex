#' @title convert to a flat data table
#' @description trims the data to an \code{IsoplotR} \code{UPb}
#'     \code{type=5} table
#' @details converts an \code{IsoplotR} \code{type=9} object to a data
#'     table that is formatted as a \code{type=5} input table
#' @param x the output of \code{calibrate}
#' @return a table n \code{IsoplotR} object of class \code{UPb}
#'     (\code{format=5})
#' @examples
#' data(Cameca,package="simplex")
#' stand <- standards(dat=Cameca,prefix='Plesovice')
#' fit <- calibration(stand=stand,oxide='UO2')
#' samp <- unknowns(dat=Cameca,prefix='Qinghu')
#' cal <- calibrate(samp,fit)
#' tab <- data2table(cal)
#' write.csv(tab,file='Qinghu.csv')
#' @export
data2table <- function(lr){
    snames <- lr$snames
    ns <- length(snames)
    geti <- function(lr,num,den){
        ns <- length(lr$snames)
        j <- which(lr$num%in%num & lr$den%in%den)
        (j-1)*ns + (1:ns)
    }
    nr <- length(lr$num)
    J <- matrix(0,nr*ns,nr*ns)
    iU <- geti(lr=lr,num='Pb206',den='U238')
    iTh <- geti(lr=lr,num='Pb208',den='Th232')
    iPb204 <- geti(lr=lr,num='Pb204',den='Pb206')
    iPb207 <- geti(lr=lr,num='Pb207',den='Pb206')
    iPb208 <- geti(lr=lr,num='Pb208',den='Pb206')
    U238Pb206 <- exp(-lr$x[iU])
    Th232Pb208 <- exp(-lr$x[iTh])
    Pb204Pb206 <- exp(lr$x[iPb204])
    Pb207Pb206 <- exp(lr$x[iPb207])
    Pb208Pb206 <- exp(lr$x[iPb208])
    J[iU,iU] <- -U238Pb206
    J[iTh,iTh] <- -Th232Pb208
    J[iPb204,iPb204] <- Pb204Pb206
    J[iPb207,iPb207] <- Pb207Pb206
    J[iPb208,iPb208] <- Pb208Pb206
    E <- J %*% lr$cov %*% t(J)
    err <- sqrt(diag(E))
    sU238Pb206 <- err[iU]
    sTh232Pb208 <- err[iTh]
    sPb204Pb206 <- err[iPb204]
    sPb207Pb206 <- err[iPb207]
    sPb208Pb206 <- err[iPb208]
    out <- cbind(Pb207Pb206,sPb207Pb206,Pb204Pb206,
                 sPb204Pb206,Pb208Pb206,sPb208Pb206)
    labels <- c('Pb207Pb206','sPb207Pb206','Pb204Pb206',
                'sPb204Pb206','Pb208Pb206','sPb208Pb206')
    if (length(U238Pb206)>0){
        out <- cbind(U238Pb206,sU238Pb206,out)
        labels <- c('U238Pb206','sU238Pb206',labels)
    }
    if (length(Th232Pb208)>0){
        out <- cbind(out,Th232Pb208,sTh232Pb208)
        labels <- c(labels,'Th232Pb208','sTh232Pb208') 
    }
    colnames(out) <- labels
    rownames(out) <- snames
    out
}
