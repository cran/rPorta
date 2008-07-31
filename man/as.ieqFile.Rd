\name{as.ieqFile}
\alias{as.ieqFile}
\title{Coerce to ieqFile object}
\description{
\code{as.ieqFile} attempts to turn its argument into an \code{ieqFile} object.
}
\usage{
as.ieqFile(v, sign=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{argument that should be turned into an \code{ieqFile} object}
  \item{sign}{vector of signs (-1 for <=, 0 for ==, and 1 for >=). Default is >=}
}
\details{
Constructs an \code{ieqFile} object consisting only of the inequalities with the help of \code{\link{as.ieq}}.}
\value{
 Returns an \code{ieqFile} object.
}
\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{
  \code{\link{as.ieq}}, \code{"\linkS4class{ieqFile}"}
}
\examples{
# x1,x2,x3>=0
as.ieqFile(cbind(diag(3),rep(0,3)))

# x1,x2,x3==0
as.ieqFile(cbind(diag(3),rep(0,3)),sign=rep(0,3))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{interface}