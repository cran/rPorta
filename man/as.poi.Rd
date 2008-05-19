\name{as.poi}
\alias{as.poi}
\title{Coerce to poi object}
\description{
\code{as.poi} attempts to turn its argument into a \code{poi} object.
}
\usage{
as.poi(v)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{argument that should be turned into a \code{poi} object}
}
\details{
Matrices with numeric values are transformed to rational values. 
This is done in a way that assumes that the numeric values in \code{v} represent rational values, e.g. \code{as.poi(0.8)=4/5}, but \code{as.poi(1/3)=333333333/1e+09}. For exact transformation, one possibility is to use the R package gmp, which gives the exact rational values, e.g. \cr \code{as.bigq(0.8)=3602879701896397/4503599627370496} \cr 
\code{as.bigq(1/3)=6004799503160661/18014398509481984} \cr as represented by the computer.

\code{as.poi} also transforms strings like "1/3" to \code{poi} objects.
}
\value{
 Returns a \code{poi} object.
}
\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{
  \code{"\linkS4class{poi}"}
}
\examples{
# Convert the numeric value 0.8 to a poi object
as.poi(0.8)

# Convert the character value "1/3" to a poi object
as.poi("1/3")

# Convert a matrix to a poi object
as.poi(matrix(1:16/4,ncol=4))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{interface}