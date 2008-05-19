\name{convertToString-method}
\alias{convertToString-method}
\alias{convertToString}
\title{Method "convertToString"}
\description{Converts some of the objects provided by rPorta to character strings, namely \code{poi}, \code{ieq}, \code{poiFile}, and \code{ieqFile}. These character strings correspond to the format used by Porta.  
}
\usage{
## S4 method for all signatures: 
convertToString(object)
}
\arguments{
  \item{object}{Object of one of the supported classes (\code{poi}, \code{ieq}, \code{poiFile}, or \code{ieqFile})}
}
\value{
Returns a character string in PORTA's format.
}
\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{
  \code{"\linkS4class{ieq}"}, \code{"\linkS4class{poi}"}, \code{"\linkS4class{ieqFile}"}, \code{"\linkS4class{poiFile}"}
}
\examples{
# Returns a String representation of a poiFile object
convertToString(example.poi())

# Returns a String representation of an ieqFile object
convertToString(example.ieq())
}
\keyword{methods}
