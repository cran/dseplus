\name{curvature}
\alias{curvature}
\alias{curvature.default}
\alias{curvature.Darray}

\title{Curvature}
\description{Curvature calculations and summary statistics as in Bates and Watts.
}
\usage{
    curvature(func, ...)
    \method{curvature}{default}(func, x, 
       method="Richardson", method.args=list(d=0.01, eps=1e-4, r=6, v=2),
        signif=0.05, show.details=FALSE, warn=TRUE, ...)
    \method{curvature}{Darray}(func, signif = 0.05,
      show.extra.details=FALSE, show.details=show.extra.details, warn=TRUE, ...)
}
\arguments{
    \item{func}{a function for the default method, a Darray object, or
       a object for which a specific method is defined (e.g. a TSestModel).}
    \item{x}{parameters to the function.}
    \item{signif}{The significance level for F test (passed to \code{qf}).}
    \item{show.details}{logical indicating if intermediate calculations 
        should be printed.}
    \item{show.extra.details}{logical indicating if extra intermediate calculations 
        should be printed.}
    \item{warn}{see \code{effectiveCurvature}.}
    \item{method}{string indicating the numerical approximation method.}
    \item{method.args}{list with arguments to \code{method} (see \code{grad}).}
    \item{...}{additional arguments passed to \code{func}.}
}
\details{
 This function is generic. It can be called with a function and parameter x,
 in which case the Bates and Watts D matrix is calculated, or with a previously
 calculated D matrix. A function \code{func} should return a sample space 
 vector at the  parameter value x.
 The method for a \code{Darray} object works on the result from genD
 which has already done most of the calculations. The Darray has an 
 element Dlist with the 3 elements as follows:
   D is a matrix of first(gradients) and second order partial
      derivatives organized in the same manner as Bates and 
      Watts. (first p columns are the gradients and the 
      next p(p-1)/2 columns are the lower triangle of the Hessian).
   p is the dimension of the parameter space=dim of the tangent space.
   f0 is the function value at the point where the matrix D 
        was calculated. (The calculation should not/does not? depend on
        this value - but it should be the right dimension and 0's do
        not work.
}
\value{
A list is returned (with invisible). Curvature summary statistics as in Bates 
and Watts are in the element stat. A representation of the Bates and 
Watts D matrix is in the element Dlist. The curvature array C as in
Batts and Watts defn. (7.16) p242 and examples p244 &  p245 is in the 
elements C.parameter and C.intrinsic.
}
\section{Side Effects}{Curvature summary statistics as in Bates and Watts are printed.}
\references{ 
Bates and Watts(1983), 'Nonlinear Regression Analysis and Its Applications.'
}
\seealso{
   \code{\link[numDeriv]{genD}}
   \code{\link{curvature.TSestModel}}
   \code{\link{effectiveCurvature}}
}
\examples{
    func <- function(x){c(x[1], x[1], x[2]^2)}
 #   curvature(func, c(2,2))
}
%\keyword{DSE}
\keyword{ts}


