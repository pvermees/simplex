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
    getij <- function(lr,num,den){
        ns <- length(lr$snames)
        i <- which(lr$num%in%num & lr$den%in%den)
        j <- (i-1)*ns + (1:ns)
        list(i=i,j=j)
    }
    nr <- length(lr$num)
    J <- matrix(0,nr,nr*ns)
    ijU <- getij(lr=lr,num='Pb206',den='U238')
    ijTh <- getij(lr=lr,num='Pb208',den='Th232')
    ijPb204 <- getij(lr=lr,num='Pb204',den='Pb206')
    ijPb207 <- getij(lr=lr,num='Pb207',den='Pb206')
    ijPb208 <- getij(lr=lr,num='Pb208',den='Pb206')
    U238Pb206 <- exp(-lr$x[ijU$j])
    Th232Pb208 <- exp(-lr$x[ijTh$j])
    Pb204Pb206 <- exp(lr$x[ijPb204$j])
    Pb207Pb206 <- exp(lr$x[ijPb207$j])
    Pb208Pb206 <- exp(lr$x[ijPb208$j])
    J[ijU$i,ijU$j] <- -U238Pb206
    J[ijTh$i,ijTh$j] <- -Th232Pb208
    J[ijPb204$i,ijPb204$j] <- Pb204Pb206
    J[ijPb207$i,ijPb207$j] <- Pb207Pb206
    J[ijPb208$i,ijPb208$j] <- Pb208Pb206
    E <- J %*% lr$cov %*% t(J)
    err <- sqrt(diag(E))
    sU238Pb206 <- err[ijU$i]
    sTh232Pb208 <- err[ijTh$i]
    sPb204Pb206 <- err[ijPb204$i]
    sPb207Pb206 <- err[ijPb207$i]
    sPb208Pb206 <- err[ijPb208$i]
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
