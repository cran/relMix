\name{relMix}
\alias{relMix}
\title{
Relationship inference based on mixtures. 
}
\description{
Calculates likelihood ratios for relationship inference involving mixtures and missing reference profiles.
}
\usage{
relMix(pedigrees, locus, datamatrix, kinship = 0)
}
\arguments{
  \item{pedigrees}{
	A list of pedigrees defined using FamiliasPedigree in Familias
}
  \item{locus}{
	A Familias locus, including mutation matrices. 
}
  \item{datamatrix}{
	Each line corresponds to one constellation of genotypes for the involved individuals 
}
  \item{kinship}{
	Defines the theta-parameter
}
}
\details{
The function requires the package Familias and calls on the function FamiliasPedigree(). 
}
\value{
The likelihood ratios for the pedigrees and detailed output for each term
}
\references{
Egeland, Kling, Mostad (2015) \pkg{Familias}
}
\author{
Navreet Kaur, Thore Egeland
}

\seealso{
Familias, see www.familias.name
}
\examples{
#Example 1
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
locus <- FamiliasLocus(frequencies=rep(0.2,5),name="V1",  
                       allelenames= c(1:5), femaleMutationRate=R,maleMutationRate=R)
relMix(pedigrees,locus,datamatrix,theta)
 
#example 2
require(Familias)
require(paramlink)
data(db)
persons <- c("CH","MO","AF")
ped1 <- pedigree(id=persons, dadid=c("AF",NA, NA), momid=c("MO", NA,NA), 
                 sex=c("male", "female", "male"))
ped2 <- pedigree(id=c(persons, "TF"), dadid=c("TF", NA, NA,NA), 
                 momid=c("MO", NA, NA,NA), sex=c("male", "female", "male", "male"))
pedigrees <- list(isFather = ped1, unrelated=ped2)
foo=db[db$Marker=="TH01",]
p=foo[c(4,5,7,8),3] 
p=c(p,1-sum(p)) #Frequencies for alleles 1,2,3,4,99 sum to 1
R=0
locus <- FamiliasLocus(frequencies=p,name="THO1.1",  
                        allelenames= c(1:4,99), femaleMutationRate=R,maleMutationRate=R)
CH=c(1,2,1,2,1,3,1,3,2,3,2,3)
MO=c(1,3,2,3,2,1,2,3,1,2,1,3)
AF=rep(c(3,4),6)
datamatrix <- data.frame(rbind(CH,MO,AF) )
rownames(datamatrix) <- c("CH","MO","AF")
relMix(pedigrees,locus,datamatrix, kinship=0.01)
}

