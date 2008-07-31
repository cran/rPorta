\name{as.ieq}
\alias{as.ieq}
\title{Coerce to ieq object}
\description{
\code{as.ieq} attempts to turn its argument into an \code{ieq} object.
}
\usage{
as.ieq(v, sign=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{argument that should be turned into an \code{ieq} object}
  \item{sign}{vector of signs (-1 for <=, 0 for ==, and 1 for >=). Default is >=}
}
\details{
Handles the same inputs as \code{\link{as.poi}} but the input has to be a matrix with at least two columns.}
\value{
 Returns an \code{ieq} object.
}
\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{
  \code{\link{as.poi}}, \code{"\linkS4class{ieq}"}
}
\examples{
# x1,x2,x3>=0
as.ieq(cbind(diag(3),rep(0,3)))

# x1,x2,x3==0
as.ieq(cbind(diag(3),rep(0,3)),sign=rep(0,3))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{interface}