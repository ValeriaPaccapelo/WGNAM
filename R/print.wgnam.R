##  Method print.wgnam
##
##' @title print method for 'wgnam'
##' @description Provides a simple printing of the selected QTL after using 'wgnam'.
##' @param object an object of class 'wgnam'
##' @param namObj nam genetic data
##' @param \ldots Additional optional objects
##' @author Valeria Paccapelo
##' @export
##'
print.wgnam <- function(object, namObj, ...)
{
  if (missing(namObj))
    stop("namObj is a required argument")
  # if (!inherits(namObj, "nam"))
  #    stop("namObj is not of class \"nam\"")
  if (is.null(object$Assoc$mark))
    cat("There are no significant putative QTL's\n")
  else {
    qtlm <- getnamAssoc(object, namObj)
    for (z in 1:nrow(qtlm)) {
      int <- paste(qtlm[z, 1], qtlm[z, 2], sep = ".")
      cat("\nPutative QTL found at marker",
               int, "\nMarker is", qtlm[z, 3], "\n")
    }
  }
}
