\name{span}
\alias{span}
\alias{span.default}

\title{Calculate Span of Tangent Plane}
\description{Calculate the dimension of the tangent space}
\usage{
    span(func, x, func.args=NULL, d=0.01, eps=1e-4, r=6,
      show.details=FALSE, ...)
    \method{span}{default}(func, x, func.args=NULL, d=0.01, eps=1e-4, r=6,
      show.details=FALSE, ...)
}
\arguments{
    \item{func}{a function which returns the residual vector for a given
         parameter vector.}
    \item{...}{arguments passed to other methods.}
    \item{x}{parameter vector first argument to function func indicating the
        point with respect to which the derivative is calculated.}
    \item{func.args}{list with additional argument to function \code{func}.}
    \item{d}{the fraction of x to use for the initial numerical approximation.
        The default is set to 0.01*x or eps if x is 0.0.}
    \item{eps}{Used instead of d for elements of x which are zero.}
    \item{r}{The number of Richardson improvement iterations (repetions with 
        successly smaller d).}
    \item{show.details}{logical indicating if detailed calculations should 
        be shown.}
}
\value{
The singular values of the matrix of tangent vectors are returned.
}
\sideffects{
If show.details is T then intermediate calculations are printed.
}
\details{
The first argument of a function must be a vector.
span performs a svd of the tangent vectors at the point x. This can be used
to calculate the dimension of the tangent space (ie. by over specifying 
the model and counting the number of significant singular values).
This function uses Richardson extrapolation  (for more details      
see the functions richardson.grad and genD) to get a numerical 
approximation of the tangent vectors to the parameter 
manifold. SVD is then used to calculate their span.
}
\seealso{
\code{\link{span.TSestModel}},
\code{\link{richardson.grad}},
\code{\link{genD}}
}
\examples{
func <- function(x){c(x[1], x[1], x[2]^2)}
span(func, c(2,2))
span(func, c(2,5))
span(func, c(2,2,5))
}
%\keyword{DSE}
\keyword{ts}

