\name{extRelMix}
\alias{extRelMix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Relationship inference based on mixtures, incorporating mutations, theta-corrections and silent alleles. 
}
\description{
The function extends the relMix-function and handles mutations, theta-correction and silent alleles for relationship inference based on mixtures.  
}
\usage{
extRelMix(E, AF, marker, R, theta, ped, silent = 0)
}
\arguments{
  \item{E}{
A DNA mixture of unknown number of contributors
}
  \item{AF}{
A known reference profile which is tested for relationship with one or more of the contributors to the DNA mixture E. 
}
  \item{marker}{
A marker found in the allele database of \code{\link{db}}
}
  \item{R}{
Mutation rate (equal for males and females)
}
  \item{theta}{
Kinship parameter (theta correction)
}
  \item{ped}{
 A list of pedigrees \code{\link{FamiliasPedigree}}
}
  \item{silent}{
Frequency of silent allele. By default 0. 
}
}
\details{
The function requires the package Familias and calls on the function FamiliasLocus() defined in Familias. 
For silent alleles; frequencies are scaled to sum to 1 if necessary. Currently only
integers can be used for allele names if silent >0
}
\value{
Returns the likelihood ratios for the pedigrees and detailed output for each term in the sum
}
\references{
Egeland, Kling, Mostad (2015) \pkg{Familias}
Kaur et al (2015)
}
\author{
Thore Egeland <thore.egeland@nmbu.no>, Navreet Kaur <navreet.kaur@nmbu.no>
}

\seealso{
Familias, see www.familias.name
}
\examples{
require(Familias)
require(paramlink)
data(db)
#Define pedigrees, only needed once
persons <- c( "CH", "MO", "AF")
ped1 <- pedigree( id = persons, dadid = c( "AF", NA, NA), 
                 momid = c( "MO", NA, NA), sex <- c( "male", "female", "male"))
ped2 <- pedigree( id = c( persons, "TF"), 
                 dadid = c( "TF", NA, NA, NA), 
                 momid =c ( "MO", NA, NA, NA), 
                 sex =c ( "male", "female", "male", "male"))
pedigrees <- list( isFather = ped1, unrelated = ped2)
res <- extRelMix(E=c(14, 16, 18), AF = c(15,18), marker = "D3S1358", R=0, theta=0, ped=pedigrees)
LR <- res$LR
LR
}