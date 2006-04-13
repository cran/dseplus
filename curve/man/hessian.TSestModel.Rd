\name{hessian.TSestModel}
\alias{hessian.TSestModel}

\title{Calculate Hessian Matrix}
\description{Calculatate the hessian matrix of a TSmodel at a parameter value.}
\usage{
    \method{hessian}{TSestModel}(func, x=coef(func), 
     method="Richardson", method.args=list(d=0.01, eps=1e-4, r=6, v=2),
     Shape=TSmodel(func), data=TSdata(func), ...)
}
\arguments{
    \item{func}{a TSestModel object which is used as a function mapping
       coefficients (parameters) to likelihood.}
    \item{x}{The parameter point at which the hessian is calculated.}
    \item{method}{string indicating the numerical approximation method.}
    \item{method.args}{list with arguments to \code{method} (see \code{grad}).}
    \item{Shape}{a TSmodel in which the parameters should be used.}
    \item{data}{TSdata to use in teh evaluation.}
    \item{...}{additional arguments passed to \code{func}.}
}
\value{a matrix of second derivative of the likelihood (Fisher Information).}
\details{
    This function calculates the second derivative of the likelihood for the
    model at its specified parameter value using given data.
}
\seealso{
    \code{\link[numDeriv]{hessian}}
    \code{\link[numDeriv]{genD}}
    \code{\link[numDeriv]{grad}}
    \code{\link{span}}
}

\concept{DSE}
\keyword{ts}

