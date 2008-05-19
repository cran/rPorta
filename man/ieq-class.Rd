\name{ieq-class}
\alias{ieq-class}
\alias{ieq}
\alias{show,ieq-method}
\alias{convertToString,ieq-method}
\alias{getNumerator,ieq-method}
\alias{getDenominator,ieq-method}
\alias{index,ieq-method}
\docType{class}
\title{Class "ieq" }
\description{Encapsulates inequalities in PORTA's ieq format.}
\section{Objects from the Class}{
An \code{ieq} object holds three slots.
}
\section{Slots}{
	 \describe{
    \item{\code{num}:}{Matrix containing the numerators of the inequalities}
    \item{\code{den}:}{Matrix containing the denominators of the inequalities}
    \item{\code{sign}:}{Vector containing either -1,0, or 1 for each inequality, corresponding to '<=', '==', and '>='}    
  }
}
\section{Methods}{
  \describe{
    \item{convertToString}{Convert to character string in PORTA's format}
    \item{getNumerator}{Returns this objects numerator matrix}
    \item{getDenominator}{Returns this objects denominator matrix}
    \item{index}{May be used to obtain subsets of the object}    
	 }
}
\section{Extends}{
Class \code{"poi"}, directly.
}
\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{
  \code{"\linkS4class{ieqFile}"}, \code{"\linkS4class{poi}"}, \code{"\linkS4class{poiFile}"}
}
\keyword{classes}
