\name{failureRegions}
\alias{failureRegions}
\alias{rConeReductionINT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Determination of unfeasible regions}
\description{
Determines unfeasible regions inside a given parameter space based on emerging failure parts,
i.e. missing values, in earlier experiments. 
}
\usage{
failureRegions(experiments, parameterspace, fail)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{experiments}{An object of class \code{poi}. The rows contain the points with known results, e.g. already solved experiments. }
  \item{parameterspace}{An object of class \code{poi} representing the parameter space as a grid. The rows should contain the points of the grid.}
  \item{fail}{A logical vector indicating which of the points of \code{experiments} are failure points. TRUE for failure points, FALSE for non-failure points.}
}
\details{
The unfeasible regions are determined as follows. Suppose the design space is
 covered roughly with a set of experiments specified in \code{experiments} and
  the results contain some failure points and some non-failure points which is indicated by \code{fail}. 
  Assuming the feasible area is convex, the space lying, from the viewpoint of the non-failure point, 
  behind the failure points has to be part of the failure region. For each non-failure point a polyhedral 
  convex cone is spanned over all non-failure points. Thereby the minimum needed polyhedral cone is constructed,
   i.e. only the extreme rays are used. Then each cone is mirrored to the opposite side of the failure points.
    Those points from \code{parameterspace} laying inside one of the mirrored cones then belong to the unfeasible region.
}
\value{
\code{failureRegions} returns an object of class \code{parameterrestriction}.  
}
\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{\code{"\linkS4class{paramspacerestriction}"}, \code{"\linkS4class{poi}"}}
\examples{
# Construct a 3-dimensional grid with values 0.0,0.2,...,1.0
parameterspaceNum<-matrix(nrow=216,ncol=3)
for (i in(0:215)) parameterspaceNum[i+1,]<-c(i\%/\%6^2,i\%/\%6\%\%6,i\%\%6)                     %/%6^2,i%/%6%%6,i%%6)
parameterspace<-new("poi",num=parameterspaceNum,den=matrix(5,nrow=216,ncol=3))

# Construct a poi object with 10 experiments
experimentsNum<-matrix(c(
	4,1,4,
	2,1,0,
	1,1,5,
	1,4,4,
	5,1,1,
	1,1,2,
	4,4,4,
	4,4,1,
	1,1,1,
	1,4,1	
),byrow=TRUE,nrow=10,ncol=3)  
experiments<-new("poi",num=experimentsNum,den=matrix(5,nrow=10,ncol=3))								

# Declare which experiments succeeded
fail <- c(0,1,0,0,0,1,0,0,1,0)

# Start
result<-failureRegions(experiments,parameterspace,fail)
show(result)

# Show summary
getParamspaceInfo(result)
}
\keyword{interface}