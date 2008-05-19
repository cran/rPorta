\name{portsort-method}
\alias{portsort-method}
\alias{portsort}
\alias{portsort,portaFile-method}
\title{Method "portsort"}
\description{
Sort inequality or point systems.
}
\usage{
## S4 method for signature 'portaFile': 
portsort(object)
}
\arguments{
  \item{object}{Object of class portaFile}
}
\details{
\code{portsort} puts the points or inequalities into an increasing order  
according to the following criteria:
\enumerate{
\item right hand sides of inequalities or equations
\item frequency  of the values -5 .. -1, 1 .. 5
\item lexicographical order
}
}
\value{
  Returns a sorted version of the input object.
}
\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{
  \code{"\linkS4class{portaFile}"}
}
\examples{
# Sorting of a poiFile object
portsort(example.poi())

# Sorting of an ieqFile object
portsort(example.ieq())
}
\keyword{methods}