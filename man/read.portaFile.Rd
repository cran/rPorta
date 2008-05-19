\name{read.portaFile}
\alias{read.portaFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parse Porta files}
\description{
\code{read.portaFile} attempts to turn the .poi or .ieq file given by its argument into a \code{portaFile} object.
}
\usage{
read.portaFile(file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{the name of the file which should be turned into a \code{portaFile} object}
}
\value{
Returns a \code{poiFile} or \code{ieqFile} object.
}
\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{
 \code{"\linkS4class{portaFile}"}, \code{"\linkS4class{poiFile}"}, \code{"\linkS4class{ieqFile}"}
}
\examples{\dontrun{	
# Could be used for one of the example files provided by PORTA
read.portaFile("example.poi")
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{interface}