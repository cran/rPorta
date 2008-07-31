\name{as.poiFile}
\alias{as.poiFile}
\title{Coerce to poiFile object}
\description{
\code{as.poiFile} attempts to turn its argument into an \code{poiFile} object.
}
\usage{
as.poiFile(v, hull=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{argument that should be turned into an \code{poiFile} object}
  \item{hull}{determines whether the convex\_hull or the convex\_cone slot should be used}
}
\details{
Constructs an \code{poiFile} object consisting only of the points of the convex cone or the convex hull with the help of \code{\link{as.poi}}.}
\value{
 Returns an \code{poiFile} object.
}
\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{
  \code{\link{as.poi}}, \code{"\linkS4class{poiFile}"}
}
\examples{
# (1,0,0), (0,1,0), and (0,0,1) as convex hull points
as.poiFile(diag(3))

# # (1,0,0), (0,1,0), and (0,0,1) as convex cone points
as.poiFile(diag(3),hull=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{interface}