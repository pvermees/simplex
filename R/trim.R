trim_pd <- function(ifname,ofname,prefixes='TEM',blocklength=137){
    fi <- file(ifname)
    file.create(ofname)
    open(fi)
    header <- readLines(fi,n=3,warn=FALSE)
    cat(paste(header,collapse='\n'),file=ofname)
    while (TRUE) {
        block <- readLines(fi,n=blocklength,warn=FALSE)
        if (length(block)<1){
            break
        } else {
            sname <- strsplit(block[[3]],split=', ')[[1]][1]
            found <- any(sapply(prefixes, grepl, sname))
            if (found){
                cat('\n',paste(block,collapse='\n'),file=ofname,append=TRUE)
            }
        }
    }
    close(fi)
}

trim_op <- function(ifname,ofname,prefixes='TEM',blocklength=38){
    fi <- file(ifname)
    file.create(ofname)
    open(fi)
    while (TRUE) {
        block <- readLines(fi,n=blocklength,warn=FALSE)
        if (length(block)<1){
            break
        } else {
            sname <- block[[1]]
            found <- any(sapply(prefixes, grepl, sname))
            if (found){
                cat(paste(block,collapse='\n'),'\n',file=ofname,append=TRUE)
            }
        }
    }
    close(fi)
}
