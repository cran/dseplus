\name{genD}
\alias{genD}
\alias{genD.default}

\title{Generate Bates and Watts D Matrix}
\description{Generate a matrix of function derivative information.}
\usage{
    genD(func, x, func.args=NULL, d=0.01, eps=1e-4, r=6)
    \method{genD}{default}(func, x, func.args=NULL, d=0.01, eps=1e-4, r=6)
}
\arguments{
    \item{func}{a function for which the first (vector) argument 
        is used as a parameter vector.}
    \item{...}{arguments to be passed to other methods.}
    \item{x}{The parameter vector first argument to \code{func}.}
    \item{d}{The fraction of x to use for the initial numerical approximation.}
    \item{eps}{Used instead of d for elements of x which are zero.}
    \item{r}{The number of Richardson improvement iterations.}
    \item{func.args}{list with additional argument to function \code{func}.}
}
\value{
A list with elements as follows:
   D is a matrix of first(gradients) and second order partial
      derivatives organized in the same manner as Bates and 
      Watts. (The first p columns are the gradients and the 
      next p(p-1)/2 columns are the lower triangle of the Hessian).
   p is the dimension of the parameter space=dim of the tangent space.
   f0 is the function value at the point where the matrix D 
        was calculated. 
   Also the arguments func, x, d, eps, r are returned in the list.
}
\details{The derivatives are calculated numerically using Richardson improvement.
 \code{genD.default} calculates a numerical approximation of the first and second
   derivatives (gradient and hessian) of \code{func} at the point x. The calculation
   is done by Richardson's extrapolation (see e.g.Linfield and Penny
   The method should be used when accuracy, as opposed to speed, is important.
  The size of d is iteratively reduced and the Richardson algorithm is applied to
    improve the accuracy of the computation.

}

\references{ 
G.R.Linfield and J.E.T.Penny,"Microcomputers in Numerical Analysis"}

\seealso{
   \code{\link{curvature}}
   \code{\link{genD.TSestModel}}
   \code{\link{gradNumerical}}
}
\examples{
    func <- function(x){c(x[1], x[1], x[2]^2)}
    z <- genD(func, c(2,2,5))
}
%\keyword{DSE}
\keyword{ts}


