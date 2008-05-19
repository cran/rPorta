\name{iespo-method}
\alias{iespo-method}
\alias{iespo}
\alias{iespo,ieqFile-method}
\title{Method "iespo"}
\description{
Enumeration of equations  and  inequalities  that
        are valid for a convex cone and a convex hull.
}
\usage{
## S4 method for signature 'ieqFile': 
iespo(object,poiObject,validity_table_out = FALSE)
}
\arguments{
  \item{object}{Object of class ieqFile}
  \item{poiObject}{Object of class poiFile}
  \item{validity_table_out}{Include a table in the  output  which indicates strong validity.}  
}
\details{
\code{iespo} is a simple enumeration routine which enumerates the subset
of equations  and  inequalities in \code{object} which are valid 
(not necessarily facet inducing) for the polyhedron given by the 
\code{poiObject}.   

}
\value{
Returns an \code{ieqFile} object.
}
\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{
  \code{"\linkS4class{ieqFile}"}, \code{"\linkS4class{poiFile}"}
}
\keyword{methods}
