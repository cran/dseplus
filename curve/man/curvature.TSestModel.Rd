\name{curvature.TSestModel}
\alias{curvature.TSestModel}
\title{Curvature for a TSestModel}
\description{Calculate curvature for a TSestModel.}
\usage{
    \method{curvature}{TSestModel}(func, x=coef(func),
     func.args=list(Shape=TSmodel(func), data=TSdata(func)),
      d=0.01, eps=1e-4, r=6, warn=TRUE, compiled=TRUE, ...)
}

\arguments{
    \item{func}{a TSestModel object which is used as a function mapping
       coefficients (parameters) to residuals.}
    \item{x}{parameter vector first argument to function func indicating the
        point with respect to which the derivative is calculated.}
    \item{func.args}{list with additional argument to function \code{func}.}
    \item{d}{The fraction of x to use for the initial numerical approximation.}
    \item{eps}{Used instead of d for elements of x which are zero.}
    \item{r}{The number of Richardson improvement iterations.}
    \item{warn}{logical. see \code{effectiveCurvature}.}
    \item{compiled}{logical, set FALSE to use curvature.default (for debugging).}
    \item{...}{arguments passed to other methods (currently ignored).}
}
\details{See the generic version of the function.}

\seealso{
    \code{\link{curvature}}
    \code{\link{effectiveCurvature}}
}
\examples{
   if(is.R()) data("eg1.DSE.data.diff", package="dse1")
   model <- estVARXls(TSdata( output=outputData(eg1.DSE.data.diff, series=1:2)),
                     max.lag=2)
   curvature(genD(model)) 
}
\concept{DSE}
\keyword{ts}

