\name{poi-class}
\alias{poi-class}
\alias{poi}
\alias{as.matrix,poi-method}
\alias{as.character,poi-method}
\alias{max,poi-method}
\alias{show,poi-method}
\alias{convertToString,poi-method}
\alias{getNumerator}
\alias{getDenominator}
\alias{getNumerator,poi-method}
\alias{getDenominator,poi-method}
\alias{index}
\alias{index,poi-method}
\docType{class}
\title{Class "poi" }
\description{Encapsulates points in PORTA's poi format.}
\section{Objects from the Class}{
A \code{poi} object holds two slots.
}
\section{Slots}{
	 \describe{
    \item{\code{num}:}{Matrix containing the numerators of the points}
    \item{\code{den}:}{Matrix containing the denominators of the points}
  }
}
\section{Methods}{
  \describe{
    \item{as.matrix}{Converts the values to a numeric matrix.}
	\item{convertToString}{Convert to character string in PORTA's format}
    \item{getNumerator}{Returns this objects numerator matrix}
    \item{getDenominator}{Returns this objects denominator matrix}
    \item{index}{May be used to obtain subsets of the object}    
		 }
}
\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{
  \code{"\linkS4class{poiFile}"}, \code{"\linkS4class{ieq}"}, \code{"\linkS4class{ieqFile}"}, \code{\link{as.poi}} 
}
\examples{
# Construct a 3-dimensional grid with values 0.0,0.2,...,1.0
parameterspaceNum<-matrix(nrow=216,ncol=3)
for (i in(0:215)) parameterspaceNum[i+1,]<-c(i\%/\%6^2,i\%/\%6\%\%6,i\%\%6)                    
parameterspace<-new("poi",num=parameterspaceNum,den=matrix(5,nrow=216,ncol=3))
}
\keyword{classes}
