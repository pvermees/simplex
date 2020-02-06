test_bg <- function(bg,tt,cc,i=1){
    p <- pars(dat$x[[i]])
    y <- log(cc)
    yp <- (bg[1] + bg[2]*tt)
    plot(rep(tt,2),c(y,yp),type='n')
    points(tt,y)
    lines(tt,yp)
}
