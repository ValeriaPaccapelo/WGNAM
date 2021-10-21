updateMPWgaim <- function(object, phenoData, mbfObj, Con.mat, attempts, ...)
{
    object$call$data <- quote(phenoData)
    if(!is.null(Con.mat)) object$call$constraints <- quote(Con.mat)
    class(object) <- "asreml"
    mbf.df <- mbfObj$mbf.df
    myby <- names(mbf.df)[1]
    myby.giv <- mbfObj$myby.giv
    object <- update(object, ...)
    fwarn <- function(object) {
        mon <- object$monitor
        prgam <- mon[4:nrow(mon), (ncol(mon) - 2)]
        prgam[abs(prgam) < .Machine$double.eps] <- NA
        pc <- abs((mon[4:nrow(mon), (ncol(mon) - 1)] - prgam)/prgam)
        ifelse(pc[!is.na(pc)] > 0.01, TRUE, FALSE)
    }
    att <- 1
    while (!object$converge) {
        object <- update(object, step = 10^{-att}, ...)
        att <- att + 1
        if (att > attempts) {
            mpwgaim:::mp.error.code("unstable")
            break
        }
    }
    att <- 1
    while (any(fwarn(object))) {
        object <- update(object, step = 10^{-att}, ...)
        att <- att + 1
        if (att > attempts) {
            mpwgaim:::mp.error.code("converge")
            break
        }
    }
    object
}
