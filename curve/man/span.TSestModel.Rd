\name{span.TSestModel}
\alias{span.TSestModel}

\title{Calculate Span of Tangent Plane}
\description{Calculate the dimension of the tangent space}
\usage{
    \method{span}{TSestModel}(func, x=coef(func),
        method="Richardson", method.args=list(d=0.01, eps=1e-4, r=6, v=2),
     show.details=FALSE, compiled=.DSEflags()$COMPILED,
     Shape=TSmodel(func), data=TSdata(func), ...)
}
\arguments{
    \item{func}{a TSestModel object which is used as a function mapping
       coefficients (parameters) to residuals.}
    \item{x}{parameter vector first argument to function func indicating the
        point with respect to which the derivative is calculated.}
    \item{show.details}{logical indicating if detailed calculations should 
        be shown.}
    \item{compiled}{use the compiled version of the code. (FALSE only for
       debugging.)}
    \item{method}{string indicating the numerical approximation method.}
    \item{method.args}{list with arguments to \code{method} (see \code{grad}).}
    \item{Shape}{a TSmodel in which the parameters should be used.}
    \item{data}{TSdata to use in teh evaluation.}
    \item{...}{additional arguments passed to \code{func}.}
}
\value{The singular values of the matrix of tangent vectors are returned.
}
\section{Side Effects}{If show.details is TRUE then intermediate calculations are printed.
}
\details{See the generic function.}
\seealso{
\code{\link{span}},
\code{\link[numDeriv]{grad}},
\code{\link[numDeriv]{genD}}
}

\concept{DSE}
\keyword{ts}
