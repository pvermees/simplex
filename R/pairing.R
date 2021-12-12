# initialise pairing between standard and sample
pairing <- function(lr,stand,type=1){
    num <- lr$method$num
    den <- lr$method$den
    stdratios <- stand$ratios
    if (type[1]==1){
        smpratios <- paste0(num,'/',den)
        matches <- match(stdratios,smpratios)
        out <- data.frame(std=stdratios,
                          smp=smpratios[matches],
                          stringsAsFactors=FALSE)
    } else {
        num_is_element <- (element(num) %in% elements())
        den_is_element <- (element(den) %in% elements())
        smp <- NULL
        versus <- NULL
        smp.c <- NULL
        if (any(!num_is_element)){ # any molecules?
            molnum <- num[!num_is_element]
            molden <- den[!num_is_element]
            nmol <- sum(!num_is_element)
            for (i in 1:nmol){
                has_mol_pair <- (den %in% molden[i])
                j <- which(num_is_element & has_mol_pair)
                smp <- c(smp,paste0(num[j],'/',den[j]))
                versus <- c(versus,paste0(molnum[i],'/',molden[i]))
                k <- which(grepl('204',num) & (den %in% num[j]))
                smp.c <- c(smp.c,paste0(num[k],'/',den[k]))
            }
        }
        if (any(!den_is_element)){
            molnum <- num[!den_is_element]
            molden <- den[!den_is_element]
            nmol <- sum(!den_is_element)
            for (i in 1:nmol){
                has_mol_pair <- (num %in% molnum[i])
                j <- which(den_is_element & has_mol_pair)
                smp <- c(smp,paste0(num[j],'/',den[j]))
                versus <- c(versus,paste0(molnum[i],'/',molden[i]))
                k <- which(grepl('204',den) & (num %in% den[j]))
                smp.c <- c(smp.c,paste0(num[k],'/',den[k]))
            }
        }
        std <- stdratios[match(smp,stdratios)]
        std.c <- stdratios[match(smp.c,stdratios)]
        out <- data.frame(std=std,smp=smp,std.c=std.c,smp.c=smp.c,
                          versus=versus,stringsAsFactors=FALSE)
        out$slope <- rep(NA,nrow(out))
    }
    out    
}
