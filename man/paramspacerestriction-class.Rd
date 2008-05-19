\name{paramspacerestriction-class}
\alias{paramspacerestriction-class}
\alias{paramspacerestriction}
\alias{show,paramspacerestriction-method}
\alias{convertToString,paramspacerestriction-method}
\alias{getUnfeasiblePoints}
\alias{getFeasiblePoints}
\alias{getParamspaceInfo}
\alias{getUnfeasiblePoints,paramspacerestriction-method}
\alias{getFeasiblePoints,paramspacerestriction-method}
\alias{getParamspaceInfo,paramspacerestriction-method}
\docType{class}
\title{Class "paramspacerestriction" }
\description{Encapsulates information returned by \code{failureRegions}.
%Warning: The points are possibly transformed in the way the algorithm is working. Since the function ??? needs truly rational values the values are possibly still transformed.
}
\section{Objects from the Class}{
An \code{paramspacerestriction} object holds five slots.
}
\section{Slots}{
	 \describe{
	\item{\code{parameterspace}:}{The original \code{poi} object representing the parameter space.}
    \item{\code{feasiblepoints}:}{Set of indices describing all points from parameterspace that do not lie inside one of the polyhedral convex cones, i.e. the feasible area.}
    \item{\code{excludingcone}:}{Set of indices describing which cone led to the exclusion of points from parameterspace}
    \item{\code{cones}:}{Set of the polyhedral convex cones constructed by \code{failureRegions}}
    \item{\code{conesources}:}{Set of indices describing the points the cones are based on}    
  }
}
\section{Methods}{
  \describe{
	\item{convertToString}{Convert to character string}
    \item{getUnfeasiblePoints}{Returns the set of unfeasible points}
    \item{getFeasiblePoints}{Returns the set of feasible points}    
    \item{getParamspaceInfo}{Returns a summary of the parameter space restriction. The first columns contain the points of the parameter space, the column named with PCC gives the set of indices in which cone the point is included, where 0 means that the point is not included in any cone. The last column PCC-Source gives the point which is the cone end.}        
	 }
}
\author{Robin Nunkesser \email{Robin.Nunkesser@tu-dortmund.de}}
\seealso{
  \code{\link{failureRegions}}
}
\keyword{classes}