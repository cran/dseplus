\name{genD.TSestModel}
\alias{genD.TSestModel}
\alias{genD.ARMA}
\alias{genD.innov}

\title{Generate Bates and Watts D Matrix}
\description{Generate a matrix of function derivative information.}
\usage{
    \method{genD}{TSestModel}(func, x=coef(func),
        func.args=list(Shape=TSmodel(func), data=TSdata(func)),
        d=0.01, eps=1e-4, r=6 )
    \method{genD}{ARMA}(func, x=coef(func),
        func.args=list(Shape=TSmodel(func), data=TSdata(func)),
        d=0.01, eps=1e-4, r=6)
    \method{genD}{innov}(func, x=coef(func),
        func.args=list(Shape=TSmodel(func), data=TSdata(func)),
        d=0.01, eps=1e-4, r=6)
}
\arguments{
    \item{func}{a TSestModel or TSmodel object which is used as a function 
       mapping coefficients (parameters) to residuals.}
    \item{x}{parameter vector first argument to function func indicating the
        point with respect to which the derivative is calculated.}
    \item{func.args}{list with additional argument to function \code{func}.}
    \item{d}{The fraction of x to use for the initial numerical approximation.}
    \item{eps}{Used instead of d for elements of x which are zero.}
    \item{r}{The number of Richardson improvement iterations.}
}
\value{
A list with three elements as follows:
   D is a matrix of first(gradients) and second order partial
      derivatives organized in the same manner as Bates and 
      Watts. (The first p columns are the gradients and the 
      next p(p-1)/2 columns are the lower triangle of the Hessian).
   p is the dimension of the parameter space=dim of the tangent space.
   f0 is the function value at the point where the matrix D 
        was calculated. 
}
\details{The derivatives are calculated numerically using Richardson improvement.}
\seealso{
\code{\link{genD}}
\code{\link{curvature}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- estVARXls(TSdata( output=outputData(eg1.DSE.data.diff, series=1:2)),
                      max.lag=2)
z <- genD(model)
}
%\keyword{DSE}
\keyword{ts}

