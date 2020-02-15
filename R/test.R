test_bg <- function(g,spot){
    p <- pars(spot)
    y <- log(cc)
    yp <- (bg[1] + bg[2]*tt)
    plot(rep(tt,2),c(y,yp),type='n')
    points(tt,y)
    lines(tt,yp)
}
