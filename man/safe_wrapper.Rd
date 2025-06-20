% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{safe_wrapper}
\alias{safe_wrapper}
\title{Safely evaluate an expression with warning on error}
\usage{
safe_wrapper(expr, msg_fmt, default = NULL)
}
\arguments{
\item{expr}{Expression to evaluate.}

\item{msg_fmt}{A format string passed to `sprintf`. The caught error message
will be inserted via `%s`.}

\item{default}{Value to return if evaluation fails. Defaults to `NULL`.}
}
\value{
The result of `expr`, or `default` on error.
}
\description{
Evaluates `expr` and returns its result. If an error occurs, a formatted
warning is issued and `default` is returned.
}
\keyword{internal}
