\name{span.TSestModel}
\alias{span.TSestModel}

\title{Calculate Span of Tangent Plane}
\description{Calculate the dimension of the tangent space}
\usage{
    \method{span}{TSestModel}(func, x=coef(func),
        func.args=list(Shape=TSmodel(func), data=TSdata(func)),
        d=0.01, eps=1e-4, r=6,
	show.details=FALSE, compiled=.DSECOMPILED)
}
\arguments{
    \item{func}{a TSestModel object which is used as a function mapping
       coefficients (parameters) to residuals.}
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
    \item{compiled}{use the compiled version of the code. (FALSE only for
       debugging.)}
}
\value{The singular values of the matrix of tangent vectors are returned.
}
\sideffects{If show.details is TRUE then intermediate calculations are printed.
}
\details{See the generic function.}
\seealso{
\code{\link{span}},
\code{\link{richardson.grad}},
\code{\link{genD}}
}

%\keyword{DSE}
\keyword{ts}
