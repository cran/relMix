\name{generate}
\alias{generate}
\title{
Generates genotypes for unknown contributors
}
\description{
Given a mixture, alleles for unknown contributors and the number of untyped contributors,
the genotypes of the unknown contributors are generated. The function is recursive.
}
\usage{
generate(R, K, x = 1)
}
\arguments{
\item{R}{
Integers representing the alleles of the mixtures
}
\item{K}{
Integers representing the alleles of the known contributors
}
 \item{x}{
The number of untyped contributors
}
}

\value{
A matrix. The number of rows is \code{x}, one row corresponds to one contributor.
The columns are the alleles, the two first for first genotype and so on.
}
\details{
Normally \code{x} is 4 or less. Computing time may be long for larger values of \code{x}.
}
\author{
Thore Egeland <Thore.Egeland@umb.no>
}

\examples{
require(paramlink)
#Given evidence R=1/2/3, known contribution K=1/2, the possible genotypes
#for 1,2 and 3 contributors are generated:
set1=generate(R=1:3,K=1:2,x=1)
set2=generate(R=1:3,K=1:2,x=2)
set3=generate(R=1:3,K=1:2,x=3)
stopifnot(all(dim(set3)==c(3,378)))
}

