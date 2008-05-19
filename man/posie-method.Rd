\name{posie-method}
\alias{posie-method}
\alias{posie}
\alias{posie,ieqFile-method}
\title{Method "posie"}
\description{
Enumeration of points that are valid for a linear
        system.
}
\usage{
## S4 method for signature 'ieqFile': 
posie(object,poiObject)
}
\arguments{
  \item{object}{Object of class ieqFile}
  \item{poiObject}{Object of class poiFile}
}
\details{
\code{posie} is a simple enumeration routine which determines the number 
of the points and direction vectors in \code{poiObject}
which are valid for the linear system in \code{object}. 
}
\value{
Returns the valid points in an \code{poiFile} object.
}
\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{
  \code{"\linkS4class{ieqFile}"}, \code{"\linkS4class{poiFile}"}
}
\keyword{methods}