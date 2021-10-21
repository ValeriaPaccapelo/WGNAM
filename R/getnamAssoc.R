getnamAssoc <- function(object, namObj)
{
    spe <- strsplit(object$Assoc$mark, split = "\\.") # mark instead of marker
    wchr <- unlist(lapply(spe, function(el) el[2]))
    wint <- as.numeric(unlist(lapply(spe, function(el) el[3])))
    mAssoc <- matrix(ncol = 4, nrow = length(wchr))
    for (i in 1:length(wchr)) {
        lhmark <- namObj$geno[[wchr[i]]]$map[wint[i]]
        mAssoc[i, 1:4] <- c(wchr[i], wint[i], names(lhmark), round(lhmark,2))
    }
    mAssoc
}
