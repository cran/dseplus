\name{hessian}
\alias{hessian}
\alias{hessian.default}

\title{Calculate Hessian Matrix}
\description{Calculatate the hessian matrix of a function at a parameter value.}
\usage{
    hessian(func, x, func.args=NULL, d=0.01, eps=1e-4, r=6)
    \method{hessian}{default}(func, x, func.args=NULL, d=0.01, eps=1e-4, r=6)
}
\arguments{
    \item{func}{ a function for which the first (vector) argument 
        is used as a parameter vector.}
    \item{x}{The parameter vector first argument to func.}
    \item{d}{The fraction of x to use for the initial numerical approximation.}
    \item{eps}{Used instead of d for elements of x which are zero.}
    \item{r}{The number of Richardson improvement iterations.}
    \item{func.args}{list with additional argument to function \code{func}.}
}
\value{The Hessian matrix of the function calculated at the point x.}
\details{
    This function uses genD and extracts the Hessian matrix..
}
\seealso{
    \code{\link{hessian.TSestModel}}
    \code{\link{genD}}
    \code{\link{span}}
}

%\keyword{DSE}
\keyword{ts}

