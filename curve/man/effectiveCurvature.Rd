\name{effectiveCurvature}
\alias{effectiveCurvature}

\title{Effective Curvature Utilitiy Used by Curvature and Project}
\description{xxx}
\usage{
    effectiveCurvature(cur, QRofD, residual, s.sqr, show.details=FALSE,
    warn=TRUE)
}
\arguments{
    \item{cur}{the relative curvature array. See \code{relCurvature}.}
    \item{QRofD}{QR decomposition of D array from Bates and Watts.}
    \item{residual}{point in sample space where the curvature should be
      calculated (possibly shifted by subtracting actual data to 
      give a residual).}
    \item{s.sqr}{sample estimate of the residual variance.}
    \item{show.details}{logical indicating if intermediate calculations 
        should be printed.}
    \item{warn}{see details.}
}
\value{x}
\details{
 Effective residual curvature from Bates and Watts p260
 Calculate the scaled RMS curvatures relative to a confidence disk 
 radius and extreme axis ratios. ref. Bates and Watts (1983) p254 eqn (7.23).
  and Bates and Watts J.R. Statist.Soc. B (1980).

 Transform the residual vector by multiply Q-transpose and sqrt(s.sqr*p).
 Calculate the p by p effective residual curvature matrix B and its 
 eigenvalues. Bates and Watts p260

 ...
   If I-B is not positive definite, where B is the effective residual 
   curvature matrix, and warn is TRUE, then a warning will indicate that the
   calculation does not seem to correspond to a local minimum.
}
\references{ 
Bates and Watts (1980), J.R. Statist.Soc. B.
Bates and Watts( 1983), 'Nonlinear Regression Analysis and Its Applications.'
}
\seealso{
    \code{\link{curvature}}
    \code{\link{project}}
    \code{\link{curvatureStats}}
    \code{\link{relCurvature}}
}
%\keyword{DSE}
\keyword{internal}


