\name{poiFile-class}
\alias{poiFile-class}
\alias{poiFile}
\alias{show,poiFile-method}
\alias{dim,poiFile-method}
\alias{convertToString,poiFile-method}
\alias{max,poiFile-method}
%\alias{writeToFile,poiFile-method}
\docType{class}
\title{Class "poiFile" }
\description{Encapsulates PORTA files describing polyhedra in PORTA's poi file format.}
\section{Objects from the Class}{
A \code{poiFile} object holds three slots.
}
\section{Slots}{
	 \describe{
    \item{\code{convex_hull}:}{poi object containing a set of points describing the polyhedron}
    \item{\code{convex_cone}:}{poi object containing a set of vectors spanning a convex cone}
    \item{\code{strong_validity}:}{Only used for return values containing a strong validity table}        
  }
}
\section{Methods}{
  \describe{
	\item{convertToString}{Convert to character string in PORTA's format}
    \item{dim}{\code{dim} computes the dimension of the convex hull and the convex cone
of a given set of points by using a gaussian elimination algorithm.
Moreover, in the case that the input system is not full dimensional, \code{dim} displays the 
equations satisfied by the system.}
%    \item{writeToFile}{...}
    
	 }
}
\section{Extends}{
Class \code{"portaFile"}, directly.
}
\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{
  \code{"\linkS4class{poi}"}, \code{"\linkS4class{ieq}"}, \code{"\linkS4class{ieqFile}"}
}
\examples{
# Show the definition of example.poi to provide an example
example.poi

# A second example
hull <-  as.poi(matrix(1:30,ncol=3))
cone <- new("poi", num=matrix(c(0,0,2), nrow=1,ncol=3), den=matrix(c(1,1,3),nrow=1,ncol=3))
poiFileExample<-new("poiFile", convex_hull=hull, convex_cone=cone)
}
\keyword{classes}