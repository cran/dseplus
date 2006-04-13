\name{curvature.TSestModel}
\alias{curvature.TSestModel}
\title{Curvature for a TSestModel}
\description{Calculate curvature for a TSestModel.}
\usage{
    \method{curvature}{TSestModel}(func, x=coef(func),
     method="Richardson", method.args=list(d=0.01, eps=1e-4, r=6, v=2),
     compiled=TRUE, warn=TRUE,
     Shape=TSmodel(func), data=TSdata(func), ...)
}

\arguments{
    \item{func}{a TSestModel object which is used as a function mapping
       coefficients (parameters) to residuals.}
    \item{x}{parameter vector first argument to function func indicating the
        point with respect to which the derivative is calculated.}
    \item{method}{string indicating the numerical approximation method.}
    \item{method.args}{list with arguments to \code{method} (see \code{grad}).}
    \item{warn}{logical. see \code{effectiveCurvature}.}
    \item{compiled}{logical, set FALSE to use curvature.default (for debugging).}
    \item{Shape}{a TSmodel in which the parameters should be used.}
    \item{data}{TSdata to use in teh evaluation.}
    \item{...}{additional arguments passed to \code{func}.}
}
\details{See the generic version of the function.}

\seealso{
    \code{\link{curvature}}
    \code{\link[numDeriv]{genD}}
    \code{\link[numDeriv]{grad}}
   \code{\link{effectiveCurvature}}
}
\examples{
   data("eg1.DSE.data.diff", package="dse1")
   model <- estVARXls(TSdata( output=outputData(eg1.DSE.data.diff, series=1:2)),
                     max.lag=2)
   curvature(genD(model)) 
}
\concept{DSE}
\keyword{ts}

