\name{example.poi}
\alias{example.poi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creates an object representing PORTA's example.poi}
\description{
  Creates an object representing one of PORTA's example files (example.poi).
}
\usage{
example.poi()
}
%- maybe also 'usage' for other objects documented here.
\value{
  Returns an \code{poiFile} object encapsulating this text file:

\preformatted{  
DIM = 3

CONV_SECTION
 3 3 0
 5/3 1 0
 1 5/2 0
CONE_SECTION
 0 0 2/3

END
}  
}
\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{\code{"\linkS4class{poiFile}"}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{interface}