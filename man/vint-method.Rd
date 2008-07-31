\name{vint-method}
\alias{vint-method}
\alias{vint}
\alias{vint,ieqFile-method}
\title{Method "vint"}
\description{
enumeration of integral inner points of  a  linear
       system
}
\usage{
## S4 method for signature 'ieqFile': 
vint(object)
}
\arguments{
  \item{object}{Object of class ieqFile}
}
\details{
\code{vint} enumerates all integral points within given bounds that  are
valid for a linear system.

In \code{object} lower and upper bounds for each component  must
be  given.  The correponding slots are \code{lower_bounds} and
\code{upper_bounds}. They have to contain exactly  dim
integers.  The  i-th  entry  of such a line gives the upper resp.
lower bound for the i-th component.
}
\value{
Returns the found integral points in a \code{poiFile} object.
}
\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{
  \code{"\linkS4class{ieqFile}"}
}
\examples{
# A very simple example
example<-as.ieqFile(cbind(diag(3),rep(0,3)))
example@lower_bounds<-c(-1,-1,-1)
example@upper_bounds<-c(1,1,1)
vint(example)
}
\keyword{methods}