\name{ieqFile-class}
\alias{ieqFile-class}
\alias{ieqFile}
\alias{show,ieqFile-method}
\alias{convertToString,ieqFile-method}
\alias{max,ieqFile-method}
%\alias{writeToFile,ieqFile-method}
\docType{class}
\title{Class "ieqFile" }
\description{Encapsulates PORTA files describing polyhedra in PORTA's ieq file format.}
\section{Objects from the Class}{
An \code{ieqFile} object holds six slots.
}
\section{Slots}{
	 \describe{
    \item{\code{valid}:}{May be used to specify a valid point for method \code{traf}}
    \item{\code{lower_bounds}:}{May be used to specify lower bounds for method \code{vint}}
    \item{\code{upper_bounds}:}{May be used to specify upper bounds for method \code{vint}}    
    \item{\code{elimination_order}:}{May be used to specify the elimination order for method \code{fmel}}
    \item{\code{inequalities}:}{\code{ieq} object describing the inequalities}    
    \item{\code{strong_validity}:}{Only used for return values containing a strong validity table}    
  }
}
\section{Methods}{
  \describe{
	\item{convertToString}{Convert to character string in PORTA's format}
%    \item{writeToFile}{...}
    
	 }
}
\section{Extends}{
Class \code{"portaFile"}, directly.
}
\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{
  \code{"\linkS4class{ieq}"}, \code{"\linkS4class{poi}"}, \code{"\linkS4class{poiFile}"}
}
\examples{
# Show the definition of example.ieq to provide an example
example.ieq

# A second example
ieqFileExample <- new("ieqFile", inequalities=new("ieq", num=matrix(1:100, ncol=10),
	den=matrix(rep(c(1,3,5,6),25), ncol=10), sign=rep(c(1,0,-1,0,1),2)))
}
\keyword{classes}