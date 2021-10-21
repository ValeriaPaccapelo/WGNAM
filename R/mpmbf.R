mpmbf <- function (phenoData, geneticData, myby, prop, res, eq.pai)
{
    int.cnt <- 2:dim(geneticData)[2]
    p <- length(int.cnt)
    ids <- as.character(geneticData[, myby])
    q <- length(ids)
    mats <- geneticData[,int.cnt]/100
    tmat <- t(mats)
    xsvd <- svd(crossprod(tmat))
    lambda <- NULL
    Psi <- NULL
    myby.giv <- NULL
    csumd <- cumsum(xsvd$d)
    cprop <- csumd/csumd[q]
    nc <- min(c(1:q)[!(cprop < prop)])
    U1 <- xsvd$u[,1:nc]
    S1 <- xsvd$d[1:nc]
    if(nc < q & res == TRUE) {
        S2 <- xsvd$d[(nc+1):q]
        if(eq.psi == TRUE) {
            trS2 <- sum(S2)
            Psi <- rep(trS2/q, q)
        }
        else {
            U2 <- xsvd$u[,(nc+1):q]
            Psi <- apply(U2, 1, function(el, S2) sum(el^2 * S2), S2=S2)
        }
        myby.giv <- data.frame(Row = 1:q, Column = 1:q, Value = 1/Psi)
        attr(myby.giv, "rowNames") <- ids
    }
    xsvd.half <- U1 %*% diag(sqrt(S1)) # t((t(U1) * sqrt(S1)))
    if(prop == 1 | res == FALSE) Gai <- U1 %*% diag(1/S1) %*% t(U1)
    else {
        Ga <- U1 %*% diag(S1) %*% t(U1) + diag(Psi)
        Gai <- solve(Ga)
    }
    xtra <- tmat  %*% Gai
    mbf.df <- as.data.frame(xsvd.half)
    names(mbf.df) <- paste("Tint.", 1:nc, sep = "")
    mbf.df <- cbind.data.frame(ids, mbf.df)
    names(mbf.df)[1] <- myby
    assign("mbf.df", mbf.df, .GlobalEnv) ## because asreml cannot find
    assign("myby", myby, .GlobalEnv)     ## these two structures
    assign("lambda", lambda, .GlobalEnv)  ## Not sure this is needed
    invisible(list(mbf.df = mbf.df, ids = ids, p = p, q = q,
        xtra = xtra, nc=nc, Psi = Psi, myby.giv = myby.giv))
}
