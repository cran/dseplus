\name{hessian.TSestModel}
\alias{hessian.TSestModel}

\title{Calculate Hessian Matrix}
\description{Calculatate the hessian matrix of a TSmodel at a parameter value.}
\usage{
    \method{hessian}{TSestModel}(func, x=coef(func),
        func.args=list(Shape=TSmodel(func), data=TSdata(func)),
        d=0.01, eps=1e-4, r=6)
}
\arguments{
    \item{func}{a TSestModel object which is used as a function mapping
       coefficients (parameters) to likelihood.}
    \item{x}{The parameter point at which the hessian is calculated.}
    \item{func.args}{list with additional argument to the TSmodel evaluation.}
    \item{d}{The fraction of x to use for the initial numerical approximation.}
    \item{eps}{Used instead of d for elements of x which are zero.}
    \item{r}{The number of Richardson improvement iterations.}
}
\value{a matrix of second derivative of the likelihood (Fisher Information).}
\details{
    This function calculates the second derivative of the likelihood for the
    model at its specified parameter value using given data.
}
\seealso{
    \code{\link{hessian}}
    \code{\link{genD}}
    \code{\link{span}}
}

%\keyword{DSE}
\keyword{ts}

