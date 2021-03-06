% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractIdbyBias.R
\name{extractIdbyBias}
\alias{extractIdbyBias}
\title{Extract ID by Bias: extracts ID value from loocv and percentile based on bias}
\usage{
extractIdbyBias(loocvObj, perc)
}
\arguments{
\item{loocvObj}{object from \code{\link{loocv_function}}.}

\item{perc}{Floating point value indicating N/100 percentile ID to look for.}
}
\value{
A numeric value of patient id number
}
\description{
Extract ID by Bias: extracts ID value from loocv and percentile based on bias
}
