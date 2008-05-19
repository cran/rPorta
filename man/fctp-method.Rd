\name{fctp-method}
\alias{fctp-method}
\alias{fctp}
\alias{fctp,ieqFile-method}
\title{Method "fctp"}
\description{
Checks inequalities for facet inducing property.
}
\usage{
## S4 method for signature 'ieqFile': 
fctp(object,poiObject)
}
\arguments{
  \item{object}{Object of class ieqFile}
  \item{poiObject}{Object of class poiFile}
}
\details{fctp performs a check whether the inequalities  given  in  \code{object} are facet inducing for the polyhedron given by
\code{poiObject}.

For all inequalities \code{fctp} does the following:

In a first step fctp checks if the inequality is valid.  If  this
is  not  the  case  \code{fctp} returns the points and rays which are not
valid.

In a second step \code{fctp} computes those valid points and rays 
which satisfy the inequality with equality and - if there are any -
returns them. 
}
\value{
  Returns one or more objects of class \code{poiFile}.
}
\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{
  \code{"\linkS4class{ieqFile}"}, \code{"\linkS4class{poiFile}"}
}
\keyword{methods}