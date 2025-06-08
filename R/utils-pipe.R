#' Null coalescing operator
#'
#' Returns `y` if `x` is `NULL`, otherwise `x`.
#'
#' @param x Value that might be `NULL`.
#' @param y Fallback value to use when `x` is `NULL`.
#' @return `x` if not `NULL`, otherwise `y`.
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

