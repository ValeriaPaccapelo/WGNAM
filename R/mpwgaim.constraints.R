mpwgaim.constraints <- function(object, phenoData, mbfObj, myby, ...)
{
    mbf.df <- mbfObj$mbf.df
    object$call$data <- quote(phenoData)
    object.sv <- update(object, start.values=TRUE, ...)
    gammas <- object.sv$gammas.table
    our.terms <- c('mbf("ints")!mbf("ints").var', paste("giv(", myby, ").giv", sep=""))
    which.terms <- gammas$Gamma %in% our.terms
    gam <- gammas[which.terms,]
    gam$con <- c(1,1)
    Con.mat <- asreml.constraints(~con, gammas=gam)
    invisible(list(object=object, Con.mat=Con.mat))
}

