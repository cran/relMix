\name{relMix-package}
\alias{relMix-package}
\docType{package}
\title{
\packageTitle{relMix}
}
\description{
\packageDescription{relMix}
}
\details{
The package contains three functions; relMix, extRelMix and generate. 
Function relMix() is the main function. 
}
\author{
\packageAuthor{relMix}.
Maintainer: \packageMaintainer{relMix}
}
\references{
Familias, www.familias.name
}
\examples{
require(Familias)
require(paramlink)
persons <- c("CH","MO","AF")
ped1 <- pedigree(id=persons, dadid=c("AF",NA, NA), momid=c("MO", NA,NA), 
                 sex=c("male", "female", "male"))
ped2 <- pedigree(id=c(persons, "TF"), dadid=c("TF", NA, NA,NA), 
                 momid=c("MO", NA, NA,NA), 
                 sex=c("male", "female", "male", "male"))
pedigrees <- list(isFather = ped1, unrelated=ped2)
E=1:3;gAF=c(3,4)
datamatrix=generate(E,K=NULL,2)
AF=rep(gAF,dim(datamatrix)[2]/2)
datamatrix=rbind(datamatrix,AF)
datamatrix=as.data.frame(datamatrix)
rownames(datamatrix)[c(1,2)] <- c("CH","MO")
R=0.00
theta=0.0
locus <- FamiliasLocus(frequencies=c(0.1,0.2,0.3,0.35,0.05),name="V1",  
                       allelenames= c(1:5), femaleMutationRate=R,maleMutationRate=R)
relMix(pedigrees,locus,datamatrix,theta)
}
