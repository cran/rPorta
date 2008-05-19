\name{fmel-method}
\alias{fmel-method}
\alias{fmel}
\alias{fmel,ieqFile-method}
\title{Method "fmel"}
\description{
Projection of linear system on subspaces xi = 0.
}
\usage{
## S4 method for signature 'ieqFile': 
fmel(object,chernikov_rule_off=FALSE,long_arithmetic=FALSE)
}
\arguments{
  \item{object}{Object of class ieqFile}
  \item{chernikov_rule_off}{Fourier-Motzkin elimination without using the rule  of Chernikov}  
  \item{long_arithmetic}{Use long integers for intermediate results.}  
}
\details{
\code{fmel} takes a system of linear inequalities as input and eliminates choosen
variables.  That is,  \code{fmel} projects the given system to the subspace
given by xi = 0, for i is contained in I, where I is the index set
of the variables that should be eliminated.

The set I and the elimination order are given in the input object in the slot \code{elimination_order}, containing exactly  dim
integers, where dim is the dimension of the problem.

A '0' as the i-th entry of \code{elimination_order} indicates that the i-th variable 
should not be eliminated, that is, i is not in I. An entry
'j', for 0 < j < dim, as the i-th entry
of \code{elimination_order} indicates that the i-th variable should be eliminated
in  the  j-th iteration. (All nonzero numbers must be different
and it must be possible to put them into an order 1,2,3,4...)
}
\value{
Returns the resulting inequalities as an \code{ieqFile} object.
}
\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{
  \code{"\linkS4class{ieqFile}"}
}
\examples{
# Fourier-Motzkin elimination of an ieqFile
fmel(example.ieq())
}
\keyword{methods}
