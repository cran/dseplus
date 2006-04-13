\name{span}
\alias{span}
\alias{span.default}

\title{Calculate Span of Tangent Plane}
\description{Calculate the dimension of the tangent space}
\usage{
    span(func, x, 
         method="Richardson", method.args=list(d=0.01, eps=1e-4, r=6, v=2),
         show.details=FALSE, ...)
    \method{span}{default}(func, x, 
         method="Richardson", method.args=list(d=0.01, eps=1e-4, r=6, v=2),
         show.details=FALSE, ...)
}
\arguments{
    \item{func}{a function which returns the residual vector for a given
         parameter vector.}
    \item{x}{parameter vector first argument to function func indicating the
        point with respect to which the derivative is calculated.}
    \item{show.details}{logical indicating if detailed calculations should 
        be shown.}
    \item{method}{string indicating the numerical approximation method.}
    \item{method.args}{list with arguments to \code{method} (see \code{grad}).}
    \item{...}{additional arguments passed to \code{func}.}
}
\value{
The singular values of the matrix of tangent vectors are returned.
}
\section{Side Effects}{
If show.details is T then intermediate calculations are printed.
}
\details{
The first argument of a function must be a vector.
span performs a svd of the tangent vectors at the point x. This can be used
to calculate the dimension of the tangent space (ie. by over specifying 
the model and counting the number of significant singular values).
This function uses Richardson extrapolation  (for more details      
see the functions gradRichardson and genD) to get a numerical 
approximation of the tangent vectors to the parameter 
manifold. SVD is then used to calculate their span.
}
\seealso{
\code{\link{span.TSestModel}},
\code{\link[numDeriv]{gradRichardson}},
\code{\link[numDeriv]{genD}}
}
\examples{
func <- function(x){c(x[1], x[1], x[2]^2)}
span(func, c(2,2))
span(func, c(2,5))
span(func, c(2,2,5))
}
%\keyword{DSE}
\keyword{ts}

