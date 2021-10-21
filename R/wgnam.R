#  Generic wgnam

#' @title Whole Genome Nested Association Mapping for one-stage QTL detection and estimation using ASReml-R.
#' @description
#' A computationally efficient whole genome approach to detecting significant QTL
#' in Nested Association Mapping Populations using the flexible one-stage linear mixed modelling functionality of ASReml-R.
#' @param baseModel a model object of class \code{asreml} usually
#' representing a phenotypic model with which to build the QTL model.
#' @param \ldots Further arguments to be passed to method
#' \code{\link{wgnam.asreml}}.

#' @return An object of class \code{wgnam} which also inherits the
#' class \code{mpwgaim} and \code{asreml} by default. The object returned is actually an
#' \code{asreml} object with the addition of
#' components from the QTL detection.
#' @export

#' @seealso \code{\link{summary.wgnam}}
#' @author Ari Verbyla
#' @references Paccapelo, M. V., Kelly, A. M., Christopher, J. T. and
#' Verbyla, A. P. (2021).  WGNAM: Whole genome nested association mapping.
#' (under review) Theoretical and Applied Genetics.
#' @references Verbyla, A. P., George, A. W., Cavanagh, C. C. and
#' Verbyla, K. L. (2014).  Whole genome QTL analysis for MAGIC.
#' Theoretical and Applied Genetics.  DOI:10.1007/s00122-014-2337-4.
#'

wgnam <- function(baseModel, ...) {
  UseMethod("wgnam")
}
