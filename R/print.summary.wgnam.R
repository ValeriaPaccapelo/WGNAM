## Function print.summary.wgnam
##
##' Print method for objects of class 'summary.wgnam'
##'
##' Prints the summary information found using 'summary.wgnam'
##' @title Print table of QTL information
##' @param x an object of class 'summary.wgnam'
##' @param \ldots other arguments passed to summary.wgnam
##' (ignored)
##' @author Ari Verbyla
##' @export
##'
print.summary.wgnam <- function(x,...) {
  print(x$summary)
  invisible()
}
