#' Method wgnam.asreml
#'
#' @title wgnam: Whole Genome Nested Association Mapping for QTL detection and estimation using ASReml-R
#' @description
#' A computationally efficient whole genome approach to detecting and estimating significant QTL
#' in Nested Association Mapping Populations using the flexible linear mixed modelling functionality of ASReml-R.
#' It includes code to handle
#' high-dimensional problems.  The size of computations in model
#' fitting is the number of genetic lines rather than the number of
#' markers * the number of founders. The use of wgnam requires a base model to
#' be fitted which includes all genetic and non-genetic terms that
#' are relevant for an analysis without marker information. The
#' marker information must be contained in an nam object
#' (\code{namObj} and these effects are added to the phenotypic data
#' (\code{phenoData}) to allow fitting of models with QTL. Note that
#' analysis using markers positions is used. All associations
#' are random effects in the final model.
#' @param baseModel a model object of class \code{asreml} usually
#' representing a base model with which to build the model.
#' @param phenoData a data frame containing the phenotypic elements
#' used to fit \code{baseModel}. This data is checked against the data
#' used in fitting the base model
#' @param namObj a list object containing the genotypic data,
#' usually an marker data object.
#' @param merge.by a character string or name of the column(s) in
#' \code{phenoData} and \code{namObj} to merge the phenotypic and
#' genotypic data sets.
#' @param aped a data frame for the inverse numerator realtionship
#' (either pedigree or marker based) of the parents in the NAM population.
#' The default is NULL.
#' @param prop specifies the proportion of the sum of the eigenvalues of the
#' genomic matrix MM'. Default is prop=1.
#' @param res is a logical, with FALSE meaning a "residual" term is
#' not included, while TRUE includes such a term. The default is FALSE,
#' @param eq.psi is a logical.  If it is TRUE a simple scaled identity residual term is
#' used otherwise a diagonal form is used.  The default FALSE.
#' @param TypeI a numerical value determining the level of
#' significance for detecting a QTL. The default is 0.05.
#' @param attempts An integer representing the number of attempts at
#' convergence for the QTL model. The default is 5.
#' @param data.name character string that represents the name of the
#' data frame for the final model fit using mpwgaim.  If no data.name
#' is specified, the file is saved in the working directory as the
#' name of the response dot data.
#' @param trace An automatic tracing facility. If \code{trace = TRUE}
#' then all \code{asreml} output is piped to the screen during the
#' analysis. If \code{trace = "file.txt"}, then output from all asreml
#' models is piped to "file.txt". Both trace mechanisms will
#' display a message if a QTL is detected.
#' @param verboseLev numerical value, either 0 or 1, determining the
#' level of tracing outputted during execution of the algorithm.  A 0
#' value will produce the standard model fitting output from the
#' fitted ASReml models involved in the forward selection. A value of
#' 1 will add a table of interval outlier statistics for each
#' iteration.
#' @param \ldots Any other extra arguments to be passed to each of
#' the \code{asreml} calls. These may also include \code{asreml.control}
#' arguments.
#'
#' @return An object of class \code{wgnam} which also inherits the
#' class \code{mpwgaim} and \code{asreml} by default. The object returned is actually an
#' \code{asreml} object with the addition of
#' components from the QTL detection.
#'
#' @seealso \code{\link{summary.wgnam}}
#'
#' @author Ari Verbyla
#'
#' @export
#'
#' @examples
#' \dontrun{
#' wgnam_qtl <- wgnam.asreml(aseModel, phenoData, namObj, merge.by = "Genotype", TypeI = 0.05)
#' }
#' @import asreml
#' @import MASS
#' @import mpwgaim

#' @references Paccapelo, M. V., Kelly, A. M., Christopher, J. T. and
#' Verbyla, A. P. (2021).  WGNAM: Whole genome nested association mapping.
#' (under review) Theoretical and Applied Genetics.
#' @references Verbyla, A. P., George, A. W., Cavanagh, C. C. and
#' Verbyla, K. L. (2014).  Whole genome QTL analysis for MAGIC.
#' Theoretical and Applied Genetics.  DOI:10.1007/s00122-014-2337-4.

wgnam.asreml <- function(baseModel, phenoData, namObj, merge.by = NULL, aped = NULL,
                         prop = 1, res = FALSE, eq.psi=FALSE, TypeI = 0.05, attempts = 5,
                         breakout = -1, data.name = NULL, trace = TRUE, verboseLev = 0, ...)
{
    if (missing(phenoData))
        stop("phenoData is a required argument.")
    if (missing(namObj))
        stop("namObj is a required argument.")
    if (!inherits(namObj, "cross"))
        stop("namObj is not of class \"cross\"")
    if (inherits(namObj, "nam"))
        nfounders <- namObj$nfounders
    if (is.null(nfounders))
        stop("nam object  does not have founder number")
    if (is.null(merge.by))
        stop("Need name of matching column to merge datasets.")
    if (is.null(gid <- namObj$pheno[, merge.by]))
        stop("Genotypic data does not contain column \"", merge.by,
            "\".")
    if (is.null(phenoData[, merge.by]))
        stop("Phenotypic data does not contain column \"", merge.by,
            "\".")
    mby <- pmatch(as.character(gid), as.character(phenoData[,
        merge.by]))
    if (all(is.na(mby)))
        stop("Names in Genotypic \"", merge.by, "\" column do not match any
names in Phenotypic \"", merge.by, "\" column.")
    if(prop == 1) {
        res <- FALSE  #### make sure no error
        eq.psi <- FALSE
    }
    whg <- gid %in% levels(phenoData[, merge.by])
    ids <- as.character(gid[whg])
    whg1 <- phenoData[, merge.by] %in% ids
    if(sum(!whg1) > 0) {
            cat("\nNote: na.method.X will be set to 'include' because there is missing
information for some terms that involve markers.
This may have consequences for other effects in your model.\n\n")
            baseModel$call$na.method.X <- "include"
    }
    mpid <- as.character(phenoData[, merge.by])
    mpid[!whg1] <- NA
    myby <- paste("mp", merge.by, sep="")
    phenoData[, myby] <- factor(mpid)
    if (is.character(trace)) {
        ftrace <- file(trace, "w")
        sink(trace, type = "output", append = FALSE)
        on.exit(sink(type = "output"))
        on.exit(close(ftrace), add = TRUE)
    }
    gdat <- lapply(namObj$geno, function(el) el$allele.data)
    if(any(lapply(gdat, function(x) dim(x)[2] %% nfounders) != 0))
        stop(" Genetic data is not consistent with the number of founders")
    nint <- lapply(gdat, function(el) 1:(ncol(el)/nfounders))
    lint <- unlist(lapply(nint, length))
    gnams <- paste("Chr", rep(names(namObj$geno), times = lint),
                   unlist(nint), sep = ".")
    eindex <- rep(1:nfounders,length(gnams))
    gnams <- paste(rep(gnams,each=nfounders), eindex, sep=".")
    geneticData <- cbind.data.frame(gid, do.call("cbind", gdat))
    names(geneticData) <- c(myby, gnams)
    whg <- geneticData[, myby] %in% levels(phenoData[, myby])
    geneticData <- geneticData[whg,]
    ids <- as.character(geneticData[, myby])
    ord.ids <- order(ids)
    ids <- ids[ord.ids]
    geneticData <- geneticData[ord.ids,]
    dnams <- names(phenoData)
    basedata <- eval(baseModel$call$data)
    bnams <- names(basedata)
    if (any(is.na(pmatch(bnams, dnams))))
        stop("Some baseModel data names do not match phenoData names")
    whn <- unlist(lapply(basedata, is.numeric))
    whn <- whn[whn][1]
    diff <- unique(abs(basedata[[names(whn)]] - phenoData[[names(whn)]]))
    if (length(diff[!is.na(diff)]) > 1)
        stop("Phenotypic data is not in same order as baseModel data.\n Try reordering
phenoData apppropriately\n")
    int.cnt <- 2:dim(geneticData)[2]
    state <- rep(1, length(int.cnt))
    names(state) <- names(geneticData[, int.cnt])
    if (!baseModel$converge) {
        cat("Warning: Base model has not converged. Updating base model\n")
        baseModel <- update(baseModel)
    }

    baseCon.mat <- NULL
    if(!is.null(baseModel$call$constraints)) baseCon.mat <- eval(baseModel$call$constraints)
    if(!is.null(baseModel$call$control$constraints)) baseCon.mat <- eval(baseModel$call$control$constraints)

    mbfObj <- nammbf(phenoData, geneticData, myby, prop, res, eq.pai, aped, nfounders)

    if(prop < 1) {
        nc <- mbfObj$nc
        nq <- length(levels(phenoData[, myby]))
        char.dr <- paste0("Dimension reduction from ", nq, " to ", nc)
        nchar.dr <- nchar(char.dr)
        eq.rep <- paste0(rep("=", nchar.dr+2), collapse="")
        cat("\n", char.dr, "\n")
        cat(eq.rep, "\n")
    }
    add.mark <- baseModel
    if(!is.null(add.mark$call$control))
        if(!is.null(add.mark$call$control$mbf)) {
            add.mark$call$control$mbf$ints <- list(key = c(myby, myby), dataFrame = "mbf.df")
        }
        else add.mark$call$control$mbf <- list(ints = list(key = c(myby, myby), dataFrame = "mbf.df"))
    else {
        if(!is.null(add.mark$call$mbf)) {
            add.mark$call$mbf$ints <- list(key = c(myby, myby), dataFrame = "mbf.df")
        }
        else add.mark$call$mbf <- list(ints = list(key = c(myby, myby), dataFrame = "mbf.df"))
    }
    int.form <- "~ idv(mbf('ints'))"
    if(prop < 1 & res == TRUE) {
        myby.giv <- mbfObj$myby.giv
        add.term <- paste(int.form, paste0("giv(", myby, ")"), sep="+")
        if(!is.null(add.mark$call$control)) {
            add.mark$call$control$ginverse[[myby]] <- myby.giv
        }
        else {
            add.mark$call$ginverse[[myby]] <- myby.giv
        }
        if(!is.null(add.mark$call$random)) add.form <- as.formula(paste(add.term, ".", sep="+"))
        if(is.null(add.mark$call$random)) add.form <- as.formula(add.term)
    }
    else {
        if(!is.null(add.mark$call$random)) add.form <- as.formula(paste0(int.form, " + ."))
        if(is.null(add.mark$call$random)) add.form <- as.formula(int.form)
    }
    Con.mat <- NULL
    if(!is.null(add.mark$call$random)) add.mark$call$random <- update.formula(add.mark$call$random, add.form)
    else add.mark$call$random <- add.form
    if(prop < 1 & res == TRUE) {
        add.con <- mpwgaim.constraints(add.mark, phenoData, mbfObj, myby)
        add.mark <- add.con$object
        Con.mat <- myCon.mat <- add.con$Con.mat
        if(!is.null(baseCon.mat)) {
            db <- dim(baseCon.mat)
            dc <- dim(Con.mat)
            dims <- dim(baseCon.mat)+dim(Con.mat)
            Con.mat <- matrix(0, dims[1], dims[2])
            dimnames(Con.mat) <- list(c(dimnames(baseCon.mat)[[1]], dimnames(Con.mat)[[1]]),
                                      c(dimnames(baseCon.mat)[[2]], dimnames(Con.mat)[[2]]))
            Con.mat[1:db[1],1:db[2]] <- baseCon.mat
            Con.mat[(db[1]+1):(db[1]+dc[1]),(db[2]+1):(db[2]+dc[2])] <- myCon.mat
            attr(Con.mat, "assign") <- c(attr(baseCon.mat, "assign"),attr(myCon.mat, "assign"))
        }
    }
    char.re1 <- "Random Effects Marker Model Iteration (1):"
    nchar.re1 <- nchar(char.re1)
    eq.rep <- paste0(rep("=", nchar.re1+2), collapse="")
    cat("\n", char.re1, "\n")
    cat(eq.rep, "\n")
    add.mark <- updateMPWgaim(add.mark, phenoData, mbfObj, Con.mat, attempts = attempts, ...)
    update <- FALSE
    which.i <- 1
    ##
    dmat <- data.frame(L0 = 0, L1 = 0, Statistic = 0, Pvalue = 0)
    vl <- cl <- oint <- blups <- list()
    mark.list <- list()
    mark <- coeff.var <- c()
    if(!is.null(aped)) {
        K.svd <- svd(aped)
        K.half <- K.svd$u %*% diag(sqrt(K.svd$d)) %*% t(K.svd$v)
        assign("K.half", K.half, .GlobalEnv)
    }
    repeat {
        if (update) {
            which.i <- which.i + 1
            add.mark$call$group[[last.mark]] <- baseModel$call$group[[last.mark]] <- marks.cnt
            mark.term <- paste0(" idv(grp(",last.mark,"))")
            if(!is.null(baseModel$call$random)) mark.form <- as.formula(paste0("~ . + ", mark.term, sep=""))
            if(is.null(baseModel$call$random)) mark.form <- as.formula(paste0("~ ", mark.term))
            mbfObj <- nammbf(phenoData, gD, myby, prop, res, eq.psi, aped, nfounders)
            mbf.df <- mbfObj$mdf.df
            myby.giv <- mbfObj$myby.giv
            Con.mat <- NULL
            baseModel$call$data <- quote(phenoData)
            char.rei <- paste0("Random Effects Marker Model Iteration (", which.i, "):")
            nchar.rei <- nchar(char.rei)
            eq.rep <- paste0(rep("=", nchar.rei+2), collapse="")
            cat("\n", char.rei, "\n")
            cat(eq.rep, "\n")
            baseModel <- updateMPWgaim(baseModel, phenoData, mbfObj, Con.mat, attempts = attempts,
                                     random. = mark.form, ...)
            add.mark$call$data <- quote(phenoData)
            if(prop < 1 & res == TRUE) {
                add.con <- mpwgaim.constraints(add.mark, phenoData, mbfObj, myby)
                add.mark <- add.con$object
                Con.mat <- myCon.mat <- add.con$Con.mat
                if(!is.null(baseCon.mat)) {
                    db <- dim(baseCon.mat)
                    dc <- dim(Con.mat)
                    dims <- dim(baseCon.mat)+dim(Con.mat)
                    Con.mat <- matrix(0, dims[1], dims[2])
                    dimnames(Con.mat) <- list(c(dimnames(baseCon.mat)[[1]], dimnames(Con.mat)[[1]]),
                                              c(dimnames(baseCon.mat)[[2]], dimnames(Con.mat)[[2]]))
                    Con.mat[1:db[1],1:db[2]] <- baseCon.mat
                    Con.mat[(db[1]+1):(db[1]+dc[1]),(db[2]+1):(db[2]+dc[2])] <- myCon.mat
                    attr(Con.mat, "assign") <- c(attr(baseCon.mat, "assign"),attr(myCon.mat, "assign"))
                }
            }
            char.rei <- paste0("Random Effects Marker plus background Model Iteration (", which.i, "):")
            nchar.rei <- nchar(char.rei)
            eq.rep <- paste0(rep("=", nchar.rei+2), collapse="")
            cat("\n", char.rei, "\n")
            cat(eq.rep, "\n")
            if(prop < 1 & res == TRUE) {
                if(!is.null(add.mark$call$control)) {
                    add.mark$call$control$ginverse[[myby]] <- myby.giv
                }
                else {
                    add.mark$call$ginverse[[myby]] <- myby.giv
                }
            }
            mark.form <- as.formula(paste0("~ . + ", mark.term, sep=""))
            add.mark <- updateMPWgaim(add.mark, phenoData, mbfObj, Con.mat, attempts = attempts,
                                   random. = mark.form, ...)
            list.coefs <- add.mark$coefficients$random
            zind <- grep("X\\.", names(list.coefs))
            list.coefs <- list.coefs[zind]
            cl[[which.i - 1]] <- list.coefs
            vl[[which.i - 1]] <- add.mark$vcoeff$random[zind]
        }
        baseLogL <- baseModel$loglik
        stat <- 2 * (add.mark$loglik - baseLogL)
        pvalue <- (1 - pchisq(stat, 1))/2
        if(stat < 0) {
            stat <- 0
            pvalue <- 1
        }
        if(pvalue < 0) pvalue <- 0
        cat("\nLikelihood Ratio Test Statistic: ", stat, ", P-value: ", pvalue,"\n\n")
        dmat[which.i, ] <- c(baseLogL, add.mark$loglik, stat, pvalue)
        if (pvalue > TypeI)
            break
        add.mark$call$data <- quote(phenoData)
        pick <- nam.pick(add.mark, phenoData, mbfObj, myby, aped, prop, res, Con.mat,
                           state, verboseLev, ...)
        state <- pick$state
        mark.list[[which.i]] <- pick$mark
        mark.x <- gsub("Chr\\.", "X.", mark.list[[which.i]])
        tmp <- strsplit(mark.x[1], split="\\.")
        mark[which.i] <- last.mark <- paste(tmp[[1]][1],tmp[[1]][2], tmp[[1]][3], sep=".")
        oint[[which.i]] <- pick$oint
        blups[[which.i]] <- pick$blups
        if(which.i==1) {
            ugtilde <- pick$ugtilde
        }
        tmp.0 <- geneticData[, mark.list[[which.i]]]
        if(!is.null(aped)) tmp.0 <- as.matrix(geneticData[, mark.list[[which.i]]]) %*% K.half
        tmp.1 <- cbind.data.frame(geneticData[, myby], tmp.0)
        names(tmp.1) <- c(myby, mark.x)
        phenoData <- cbind.data.frame(ord = 1:nrow(phenoData),
                                      phenoData)
        phenoData <- merge(phenoData, tmp.1, by = myby,
                           all.x = TRUE, all.y = FALSE)
        phenoData <- phenoData[order(phenoData$ord), ]
        phenoData <- phenoData[, -2]
        gbind <- (2:dim(geneticData)[2])[!as.logical(state)]
        gD <- geneticData[, -gbind]
        marks.cnt <- grep(paste(last.mark, ".", sep=''), names(phenoData), fixed=TRUE)
        if(breakout > 0) {
            if(breakout == which.i) {
                break}
        }
        update <- TRUE
    }
    sigma2 <- add.mark$sigma2
    if (names(add.mark$gammas.con[length(add.mark$gammas.con)]) ==
        "Fixed")
        sigma2 <- 1
    add.mark$Assoc$founders <- namObj$founders
    add.mark$Assoc$nfounders <- nfounders
    add.mark$Assoc$diag$dmat <- dmat
    add.mark$Assoc$diag$ugtilde <- ugtilde
    add.mark$Assoc$diag$blups <- blups
    add.mark$Assoc$diag$oint <- oint
    add.mark$Assoc$aped <- aped
    add.mark$Assoc$prop <- prop
    add.mark$Assoc$res <- res
    add.mark$Assoc$eq.psi <- eq.psi
    if (length(mark)) {
        add.mark$Assoc$mark <- mark
        add.mark$Assoc$mark.list <- mark.list
        if(breakout != 1) {
            add.mark$Assoc$effects <- cl[[which.i-1]]
            add.mark$Assoc$veffects <- vl[[which.i-1]]
            add.mark$Assoc$diag$vl <- vl
            add.mark$Assoc$diag$cl <- cl
        }
        add.mark$Assoc$state <- state
        add.mark$Assoc$Con.mat <- Con.mat
        add.mark$Assoc$mbfObj <- mbfObj
        add.mark$Assoc$myby <- myby
    }
    add.mark <- mpwgaim:::mp.envDestruct(add.mark, keep = c("phenoData", "mbfObj", "Con.mat",
        "ftrace", "data.name"))
    class(add.mark) <- c("wgnam", "wgaim", "asreml")
    if(is.null(data.name)) {
        data.name <- paste(as.character(add.mark$call$fixed[2]),
                           "data", sep = ".")
    }
    add.mark$call$data <- as.name(data.name)
    cat(" Saving final data frame as ", data.name, " in working directory\n")
    assign(data.name, phenoData, envir = parent.frame())
    add.mark
}

