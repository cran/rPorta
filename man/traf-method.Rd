\name{traf-method}
\alias{traf-method}
\alias{traf}
\alias{traf,portaFile-method}
\title{Method "traf"}
\description{
Transformation of polyhedron representations.
}
\usage{
## S4 method for signature 'portaFile': 
traf(object, opt_elim = FALSE, chernikov_rule_off = FALSE, 
	validity_table_out = FALSE,long_arithmetic=FALSE)
}
\arguments{
  \item{object}{Object of class portaFile}
  \item{opt_elim}{Use a heuristic to eliminate that variable  next,  for which the number of new inequalities is minimal (local 
criterion). 
%If this option  is set, inequalities  which are  recognized  to  be facet-inducing  for the finite 
%linear system are printed into a  file as soon as they are identified.
}
  \item{chernikov_rule_off}{Fourier-Motzkin elimination without using the rule  of Chernikov}
  \item{validity_table_out}{Include a   table in the  output  which indicates strong validity.}
  \item{long_arithmetic}{Use long integers for intermediate results.}
}
\details{
\code{traf} transforms polyhedra between the representations \code{poiFile} (convex hull of points + convex cone of vectors) and \code{ieqFile} (system of linear equations and inequalities). 
The direction  of  transformation  is  determined  by  the  class of \code{object}. All computations are carried out in rational arithmetic  to  have
guaranteed  correct  numerical  results. Rational arithmetic uses only integer operations. 
% A possible arithmetic overflow is recognized. In
%this case the computations can be restarted with a special arithmetic 
%allowing the integers to have arbitrary length. This 
%arithmetic is not as efficient as the system's integer arithmetic with
%respect to time and storage requirements.

The computation of  the  ieq-representation  is  performed  using
Gaussian  and Fourier-Motzkin elimination. In the output inequalities the
right hand sides are 0, or  determined  by  the  smallest  integer
value  for which the coefficients of the inequality are integral.
If this is not possible with system integer arithmetic 
%or if multiple precision  integer arithmetic is set
, the right hand sides are 0 or 1 or -1 and the values are reduced as far as possible.
The resulting inequalities are all facet-defining for your polyhedron and give
together with equations a minimal linear description
of your polyhedron.


If an 'ieq'-representation is given as input  and  if  0  is  not
valid for the linear system, \code{traf} needs a valid point. Such a valid point may be specified in the slot
\code{valid} of objects of class \code{ieqFile}. \code{traf} transforms  the
ieq  representation  to the poi-representation, after elimination
of equations and 0-centering, by applying the 'poi'-to-'ieq'
direction to the polar polyhedron. 

\bold{Hint}: If you give a valid
point  or if 0 is valid, then this vector may appear again in the 
resulting system, even if this vector might be redundant in a 
minimal description. (All other vectors are non-redundant.)
}
\value{
  Returns an object of class \code{poiFile} or \code{ieqFile} depending on the direction of transformation.
}
\references{K. Fukuda and A. Prodon (1996). Double Description Method Revisited. In: \emph{Combinatorics and Computer Science.} Vol. 1120 of Lecture Notes in Computer Science. Springer, London, 91--111.
}

\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{
  \code{"\linkS4class{portaFile}"}, \code{"\linkS4class{poiFile}"}, \code{"\linkS4class{ieqFile}"} 
}
\examples{
# Convert a poiFile object to an ieqFile object
traf(example.poi())

# Convert an ieqFile object to a poiFile object
traf(example.ieq())
}
\keyword{methods}