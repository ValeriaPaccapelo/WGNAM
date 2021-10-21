mat2sparse <- function (X, rowNames = dimnames(X)[[1]])
{
    which <- (X != 0 & lower.tri(X, diag = TRUE))
    df <- data.frame(row = t(row(X))[t(which)], column = t(col(X))[t(which)],
        value = t(X)[t(which)])
    if (is.null(rowNames))
        rowNames <- as.character(1:nrow(X))
    attr(df, "rowNames") <- rowNames
    df
}


