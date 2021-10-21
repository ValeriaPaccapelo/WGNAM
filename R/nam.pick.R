nam.pick <-
function (asr, phenoData, mbfObj, myby, aped, prop, res, Con.mat, state, verboseLev, ...)
{
    asr$call$data <- quote(phenoData)
    if(!is.null(Con.mat)) {
        asr$call$constraints <- Con.mat
    }
    mbf.df <- mbfObj$mbf.df
    Psi <- mbfObj$Psi
    ids <- mbfObj$ids
    xtra <- mbfObj$xtra  ###### A large matrix in one direction
    q <- mbfObj$q
    sigma2 <- asr$sigma2
    if (names(asr$gammas.con[length(asr$gammas.con)]) == "Fixed")
        sigma2 <- 1
    only.terms <- "mbf(\"ints\")"
    if(prop < 1 & res == TRUE) {
        only.terms <- c(only.terms, paste0("giv(", myby,")"))
    }
    cat(" Predict step for outlier statistics \n")
    cat("=====================================\n")
    pv <- predict(asr, classify = myby, only = only.terms, levels = ids,
                  vcov = TRUE, maxiter = 1, data=phenoData, ...)
    gamma <- pv$gammas[grep("ints", names(asr$gammas))]
    ugtilde <- pv$predictions$pvals[, "predicted.value"]
    pev <- pv$predictions$vcov
    tmp <- mbf.df[,-1]
    mbf.mat <- as.matrix(tmp)
    var.u <- sigma2 * gamma * (mbf.mat %*% t(mbf.mat))
    if(prop < 1 & res == TRUE) {
        var.u <- var.u + sigma2 * gamma * diag(Psi)
    }
    atilde <- xtra %*% ugtilde
    vuatilde <- var.u - pev
    gnams <- names(state)[as.logical(state)]
    names(atilde) <- gnams
    split.gnams <- strsplit(gnams, split = "\\.")
    mgnams <- unique(unlist(lapply(split.gnams, function(x) paste(x[1],
        x[2], x[3], sep = "."))))
    mgnams1 <- paste0(mgnams, ".")
    nisumtj2 <- c()
    noint <- length(mgnams)
    cat("Before outlier statistics\n")
    nisumtj2 <- lapply(mgnams1, function(el, atilde, vuatilde, gnams, aped, xtra) {
        wh.m <- grep(el, gnams, fixed = TRUE)
        vmat <- xtra[wh.m, ] %*% vuatilde %*% t(xtra[wh.m,])
        dvatilde <- diag(vmat)
        if(is.null(aped)) t2 <- sum(atilde[wh.m]^2)/sum(diag(vmat))
        else {
            apedi <- solve(aped)
            num <- sum(atilde[wh.m] * t(atilde[wh.m] * apedi))
            denom <- sum(apedi * vmat)
            t2 <- num/denom
        }
        list(t2=t2, dvatilde=dvatilde)
    }, atilde=atilde, vuatilde=vuatilde, gnams=gnams, aped=aped, xtra=xtra)
    cat("After outlier statistics\n")
    t2 <- unlist(lapply(nisumtj2, function(el) el$t2))
    names(t2) <- mgnams
    t2[is.na(t2)] <- 0
    markint <- names(t2)[t2 == max(t2)]
    print(markint)
    qsp <- unlist(strsplit(markint, split = "\\."))
    wint <- as.numeric(qsp[3])
    print(wint)
    wchr <- qsp[2]
    print(wchr)
    oint <- c(t2)
    blups <- state
    mark.interval <- paste("Chr", wchr, wint, sep = ".")
    mark.interval <- paste(mark.interval, ".", sep = "")
    mark <- nint <- gnams[grep(mark.interval, gnams, fixed = TRUE)]
    tint <- state[unique(unlist(lapply(split.gnams, function(x) paste(x[1],
        x[2], x[3], 1, sep = "."))))]
    tint[as.logical(tint)] <- oint
    oint <- tint
    dvatilde <- unlist(lapply(nisumtj2, function(el) el$dvatilde))
    blups[as.logical(state)] <- atilde/sqrt(dvatilde)
    if (verboseLev > 0) {
        cgen <- "Marker"
        cat(cgen, "outlier statistics \n")
        cat("=============================================== \n")
        for (i in 1:length(t2)) cat(cgen, names(t2)[i],
            "Outlier Statistic ", t2[i], "\n")
        cat("=============================================== \n\n")
    }
    state[nint] <- 0
    message("\n Found Association on chromosome ", wchr, " at marker ",
        wint, "\n")
    list(state = state, mark = mark, oint = oint, blups = blups, ugtilde=ugtilde)
}
