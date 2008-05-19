\name{portaFile-class}
\alias{portaFile-class}
\alias{show,portaFile-method}
\alias{convertToString,portaFile-method}
\alias{writeToFile,portaFile-method}
\alias{writeToFile}
\alias{rPortaInterface}
\docType{class}
\title{Class "portaFile" }
\description{Encapsulates PORTA files in PORTA's ieq and poi file format. Is the abstract superclass of \code{"ieqFile"} and \code{"poiFile"}.}
\section{Methods}{
  \describe{
	\item{convertToString}{Convert to character string in PORTA's format}
    \item{writeToFile}{Writes the character string to a file}
	 }
}
\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{
  \code{"\linkS4class{ieqFile}"}, \code{"\linkS4class{poiFile}"}
}
\keyword{classes}
