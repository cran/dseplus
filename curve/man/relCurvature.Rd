\name{relCurvature}
\alias{relCurvature}

\title{Relative Curvature Utility Used by Curvature and Project}
\description{xxx}
\usage{
    relCurvature(s.sqr, R11, R2, show.extra.details=FALSE,
          eps=sqrt(.Machine$double.eps))
}
\arguments{
    \item{s.sqr}{sample estimate of the residual variance.}
    \item{R11}{submatrix of R. See details.}
    \item{R2}{submatrix of R. See details.}
    \item{show.extra.details}{logical indicating if extra intermediate calculations 
        should be printed.}
    \item{eps}{singular values smaller than eps times the largest singular 
        value are considered to be zero.}
}
\value{relative curvature array}
\details{
 The result is the relative curvature array Bates & Watts p242-244 eqn. (7.16)
  and examples p244-245.
  
 R11 is a p by p sub matrix of R. See Bates & Watts p236
  and  R2 is the m by pp sub matrix  R12 
				     R22   
...
   If I-B is not positive definite, where B is the effective residual 
   curvature matrix, and warn is TRUE, then a warning will indicate that the
   calculation does not seem to correspond to a local minimum.
}
\seealso{
    \code{\link{curvature}}
    \code{\link{project}}
    \code{\link{effectiveCurvature}}
    \code{\link{curvatureStats}}
}
%\keyword{DSE}
\keyword{internal}



