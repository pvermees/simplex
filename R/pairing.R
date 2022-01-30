#' @title Pair isotopic ratios for multi-element calibration
#' @description Pairs oxide/element ratios with element/element ratios
#'     for multi-element calibration, typically in the context of U-Pb
#'     and Th-Pb geochronology
#' @param lr an object of class \code{logratios}
#' @param stand an object of class \code{stand}
#' @return an object of class \code{pairing}, i.e. a data frame with
#'     the following columns:
#'
#' \code{X}: the parent oxide/parent ratio (e.g., \code{'UO/U238'},
#' \code{'ThO/Th'}),
#'
#' \code{Y}: the parent (e.g. \code{'U238'}, \code{'Th232'}),
#'
#' \code{C} (optional): the common Pb ratio
#'     (e.g. \code{'Pb204/Pb206'}, \code{'Pb204/Pb208'})
#'
#' \code{slope}: either \code{'auto'} or a numerical value if the
#' calibration is to use a fixed slope.
#' @examples
#' \dontrun{
#' data('SHRIMP_UPb',package='simplex')
#' st <- standard(preset='Temora')
#' dc <- drift(x=SHRIMP_UPb)
#' lr <- logratios(x=dc)
#' cal <- calibration(lr=lr,stand=st,prefix="TEM")
#' plot(cal)
#' }
#' @export
pairing <- function(lr,stand){
    if (missing(stand)) stand <- skeletonstand(lr)
    num <- lr$method$num
    den <- lr$method$den
    OP <- NULL
    DP <- NULL
    if (!stand$measured) CD <- NULL
    num_is_element <- (element(num) %in% elements())
    den_is_element <- (element(den) %in% elements())
    if (all(num_is_element) & all(den_is_element)){
        OP <- ''
        DP <- ''
        CD <- ''
    } else {
        if (any(!num_is_element)){ # any molecules?
            molnum <- num[!num_is_element]
            molden <- den[!num_is_element]
            nmol <- sum(!num_is_element)
            for (i in 1:nmol){
                has_mol_pair <- (den %in% molden[i])
                j <- which(num_is_element & has_mol_pair)
                DP <- c(DP,paste0(num[j],'/',den[j]))
                OP <- c(OP,paste0(molnum[i],'/',molden[i]))
                if (!stand$measured){
                    k <- which(grepl('204',num) & (den %in% num[j]))
                    CD <- c(CD,paste0(num[k],'/',den[k]))
                }
            }
        }
        if (any(!den_is_element)){
            molnum <- num[!den_is_element]
            molden <- den[!den_is_element]
            nmol <- sum(!den_is_element)
            for (i in 1:nmol){
                has_mol_pair <- (num %in% molnum[i])
                j <- which(den_is_element & has_mol_pair)
                DP <- c(DP,paste0(num[j],'/',den[j]))
                OP <- c(OP,paste0(molnum[i],'/',molden[i]))
                if (!stand$measured){
                    k <- which(grepl('204',den) & (num %in% den[j]))
                    CD <- c(CD,paste0(num[k],'/',den[k]))
                }
            }
        }
    }
    out <- data.frame(X=OP,Y=DP,stringsAsFactors=FALSE)
    if (!stand$measured) out$C <- CD
    out$slope <- rep('auto',nrow(out))
    class(out) <- append(class(out),'pairing')
    out
}
