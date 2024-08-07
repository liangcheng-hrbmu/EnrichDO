% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EnrichDO_visualization.R
\name{writeDoTerms}
\alias{writeDoTerms}
\title{writeDoTerms}
\usage{
writeDoTerms(doterms = doterms, file = "doterms.txt")
}
\arguments{
\item{doterms}{a data frame of do terms.}

\item{file}{the address and name of the output file.}
}
\value{
text
}
\description{
Output DOterms as text
}
\examples{
demo.data=c(1636,351,102,2932,3077,348,4137,54209,5663,5328,23621,3416,3553)
sample3<-doEnrich(interestGenes=demo.data)
writeDoTerms(sample3)
}
\author{
Haixiu Yang
}
