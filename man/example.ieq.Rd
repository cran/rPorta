\name{example.ieq}
\alias{example.ieq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creates an object representing PORTA's example.ieq}
\description{
  Creates an object representing one of PORTA's example files (example.ieq).
}
\usage{
example.ieq()
}
%- maybe also 'usage' for other objects documented here.
\value{
  Returns an \code{ieqFile} object encapsulating this text file:

\preformatted{  
DIM = 5
VALID
3 3 0 2 3

LOWER_BOUNDS
0 1 2 2 2
UPPER_BOUNDS
2 2 2 5 5

ELIMINATION_ORDER
2 0 1 0 3

INEQUALITIES_SECTION
(  1) +27x1-28x2+57x4-37x5 == 0
(  2)           -  x4+  x5 == 1

(  1) + x2   -2x5 <=  -3
(  2)     -x3     <=   0
-2x2   + x5 <=   0
-4/15x2-1/15x5 >=  -1

END
}  
}
\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{\code{"\linkS4class{ieqFile}"}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{interface}