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
    out
}

skeletonpairing <- function(){
    
}
