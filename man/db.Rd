% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db.R
\docType{data}
\name{db}
\alias{db}
\title{Allele database}
\format{
A data frame with 324 observations on the following 3 variables:
 \describe{
\item{\code{Marker}}{a factor with levels corresponding to name of markers}
\item{\code{Allel}}{a numeric vector denoting allele}
\item{\code{Frequency}}{a numeric vector in (0,1)}}
}
\source{
Dupuy et al. (2013), unpublished.
}
\usage{
data(db)
}
\description{
Norwegian database with 17 EXS17 markers and 6 additional markers.
}
\examples{
data(db)
#Checks that frequencies add to 1
lapply(split(db$Frequency,db$Marker),sum)
#Finds number of alleles for all markers
unlist(lapply(split(db$Frequency,db$Marker),length))
#A closer look at the marker SE33
SE33=db[db$Marker=="SE33",]
barplot(SE33$Frequency)
}
\keyword{datasets}
