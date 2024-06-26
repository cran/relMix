% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkReferenceFile.R
\name{checkReferenceFile}
\alias{checkReferenceFile}
\title{Check a reference profiles file}
\usage{
checkReferenceFile(filename, mix)
}
\arguments{
\item{filename}{Path of the reference profiles file}

\item{mix}{Data frame with mixture data}
}
\value{
A list containing
\itemize{
\item \code{df} The loaded data frame, NULL if errors are present
\item \code{warning} A list of strings describing the errors that ocurred but could be fixed or that do not prevent the execution of the program.
\item \code{error} A list of strings describing the errors that ocurred that made it impossible to return a valid data frame.
If this list is not empty, then the data frame item will be null.}
}
\description{
Given a reference profile file name the function attempts to load it and compare it to the mixture file to detect possible errors.
}
\details{
See the relMix vignette for a description of the format of the reference file. The data frame with mixture data is used to compare
If warnings are found, the function attempts to fix them and explains what it has done in the warning messages.
If an error is found, checking stops and a NULL dataframe is returned. The error is described in the error messages.
}
\examples{
#Load a mixture file
mixfile <- system.file("extdata","mixture.txt",package="relMix")
mix <- checkMixtureFile(mixfile);
#Note: the mixture dataframe is passed as an argument. If the previous check failed,
#the program should not continue with the reference file check
reffile <- system.file("extdata","references.txt",package="relMix")
checkReferenceFile(reffile, mix$df);
}
\seealso{
\code{\link{checkMixtureFile}} for information on how to load a mixture file.
}
\author{
Elias Hernandis
}
