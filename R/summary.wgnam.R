#  Method summary.wgnam
#
#' @title Summary of the QTL analysis after using wgnam
#'
#' @description Performs a summary after a QTL analysis using \code{wgnam}. QTL effects are organised in a presentable format.
#'
#' @param object an object of class \code{wgnam}
#' @param namObj a list object containing the genotypic data,
#' @param merge.by The genetic line variable used in merging
#' phenotypic and genetic data in the \code{wgnam} analysis.
#' @param \ldots Optional additional objects passed to or from other
#' methods
#' @return A object of class \code{summary.wgnam}, consisting of a list
#' with a data frame of QTL information.
#' @author Ari Verbyla
#' @export
#' @examples
#' \dontrun{
#' SumNAM <- summary(wgnam_qtl)
#' }
#'
summary.wgnam <-
function (object, namObj, merge.by = "id", ...)
{
    phenoData <- eval(object$call$data)
    if (missing(namObj))
        stop("namObj is a required argument")
    # if (!inherits(namObj, "nam"))
    #     stop("namObj is not of class \"nam\"")
    if (is.null(markers <- object$Assoc$mark)) {
        cat("There are no significant marker associations\n")
        return()
    }
    else {
        sigma2 <- object$sigma2
        if (names(object$gammas.con[length(object$gammas.con)]) ==
            "Fixed")
            sigma2 <- 1
        nfounders <- object$Assoc$nfounders
        Con.mat <- object$Assoc$Con.mat
        if(!is.null(Con.mat)) object$call$constraints <- Con.mat
        assign("Con.mat", Con.mat, .GlobalEnv)
        if(!is.null(object$Assoc$aped)) {
            aped <- object$Assoc$aped
            K.svd <- svd(aped)
            K.half <- K.svd$u %*% diag(sqrt(K.svd$d)) %*% t(K.svd$v)
            assign("K.half", K.half, .GlobalEnv)
        }
        LOGP <- LogProb <- c()
        Pvalue <- c()
        effects <- i.LOGP <- i.LogProb <- list()
        i.Pvalue <- list()
        for (ii in 1:length(markers)) {
            marker <- object$Assoc$mark[ii]
            whq <- unlist(eval(object$call$group[[marker]]))
            n.int <- length(whq)
            gmarker <- gsub("X.", "", marker)
            char.dr <- paste0(" Predict step for Association of marker ", gmarker)
            nchar.dr <- nchar(char.dr)
            eq.rep <- paste0(rep("=", nchar.dr+2), collapse="")
            cat("\n", char.dr, "\n")
            cat(eq.rep, "\n")
            pred.term <- paste0("grp(", marker, ")")
            qlev <- diag(n.int)
            qlist <- list(as.vector(qlev))
            names(qlist) <- marker
            pred <- predict(object, classify = pred.term, only = paste0("grp(", marker, ")"),
                            levels = qlist, vcov = TRUE, maxiter = 1)$predictions
            effects[[ii]] <- pred$pvals$predicted.value
            vcov <- pred$vcov
            if(!is.null(object$Assoc$aped)) {
                effects[[ii]] <- as.vector(K.half %*% effects[[ii]])
                vcov <- K.half %*% vcov %*% K.half
            }
            d2.ii <- sum(effects[[ii]] %*% MASS::ginv(vcov) * effects[[ii]])
            i.d2.ii <- (effects[[ii]]^2)/diag(vcov)
            Pvalue[ii] <- pchisq(d2.ii, df = nfounders - 1, lower.tail = FALSE)
            LOGP[ii] <- round(0.5 * log(exp(d2.ii), base = 10),
                2)
            LogProb[ii] <- round(-log(Pvalue[ii], base = 10),
                2)
            i.Pvalue[[ii]] <- pchisq(i.d2.ii, df = 1, lower.tail = FALSE)/2
            i.LOGP[[ii]] <- round(0.5 * log(exp(i.d2.ii), base = 10),
                3)
            i.LogProb[[ii]] <- round(-log(i.Pvalue[[ii]], base = 10),
                2)
        }
        iPvalues <- unlist(i.Pvalue)
        iLOGPs <- unlist(i.LOGP)
        iLogProbs <- unlist(i.LogProb)

        markermat <- getnamAssoc(object, namObj)

        Pvalues <- addons <- LOGPs <- LogProbs <- rep("", length(object$Assoc$effects))
        seqs <- seq(from = 1, to = length(Pvalues), by = nfounders)
        Pvalues[seqs] <- round(Pvalue, 3)
        LOGPs[seqs] <- round(LOGP, 2)
        LogProbs[seqs] <- round(LogProb, 2)
        markermat.expand <- matrix("", nrow = dim(markermat)[1] * nfounders,
            ncol = dim(markermat)[2])
        markermat.expand[seqs, ] <- markermat
        wchr <- rep(markermat[, 1], each = nfounders)
        wint <- rep(as.numeric(markermat[, 2]), each = nfounders)
        gdat <- lapply(namObj$geno, function(el) el$allele.data)
        if(any(lapply(gdat, function(x) dim(x)[2] %% nfounders) != 0))
        stop(" Genetic data is not consistent with the number of founders")
        nint <- lapply(gdat, function(el) 1:(ncol(el)/nfounders))
        lint <- unlist(lapply(nint, length))
        gnams <- paste("Chr", rep(names(namObj$geno), times = lint),
                   unlist(nint), sep = ".")
        eindex <- rep(1:nfounders,length(gnams))
        gnams <- paste(rep(gnams,each=nfounders), eindex, sep=".")
        geneticData <- do.call("cbind", gdat)
        colnames(geneticData) <- gnams                                        # colnames instead of names
        gD <- geneticData[, !(gnams %in% unlist(object$Assoc$mark.list))]
        gnam <- unlist(lapply(strsplit(names(object$gammas),
            "!"), function(el) el[1]))
        gnam <- names(namObj$pheno)[names(namObj$pheno) %in%
            gnam]
        marker.var <- c()
        for (ii in 1:length(object$Assoc$mark.list)) {
            marks <- geneticData[, object$Assoc$mark.list[[ii]]]
            av.mark <- apply(marks, 2, mean)
            mark.mat <- av.mark %*% t(av.mark)
            marker.var[ii] <- sum((effects[[ii]] %*% (diag(av.mark) -
                mark.mat)) * effects[[ii]])
        }
        coeff.other <- sum(apply(gD, 2, mean)^2)
        other.var <- sigma2 * object$gammas[grep(paste("!", gnam,
            sep = ""), names(object$gammas))]
        if (length(other.var) == 0)
            other.var <- sigma2
        coeff.other <- c(coeff.other, rep(1, length(other.var)))
        other.var <- c(sigma2 * object$gammas[grep("ints", names(object$gammas))]/10000,
            other.var)
        nonmarker.var <- coeff.other * other.var
        full.var <- sum(marker.var) + sum(nonmarker.var)
        perc.var <- 100 * marker.var/full.var
        names(perc.var) <- markers
        adds <- round(perc.var, 1)
        addons[seqs] <- adds
        effects <- unlist(effects)
        markermat <- cbind(markermat.expand, rep(1:nfounders, length(markers)),
            round(effects, 3), round(iPvalues, 3), round(iLogProbs,
                2), Pvalues, addons, LogProbs)
        collab <- c("Chromosome", "Marker No.", "Marker", "dist (cM)",
                    "Founder", "Size", "Allele Prob",
                    "Allele LOGP", "Prob", "% var", "LOGP")
        Allelepos <- 5
        markermat <- markermat[order(wchr, wint, as.numeric(markermat[,
                                                                Allelepos])), ]
        dimnames(markermat) <- list(as.character(1:length(effects)),
            collab)
        Founder <- rep(namObj$founders, length(markers))
        markermat[, "Founder"] <- Founder
        markermat <- markermat[, -2]
        markermat <- as.data.frame(markermat)
        rownames(markermat) <- 1:dim(markermat)[1]
    }
    outObject <- list(summary = markermat,
        markers = markers, effects=effects, Prob=Pvalue, LOGP = LOGP)
    class(outObject) <- "summary.wgnam"
    invisible(outObject)
}


