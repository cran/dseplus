\name{span.TSestModel}
\alias{span.TSestModel}

\title{Calculate Span of Tangent Plane}
\description{Calculate the dimension of the tangent space}
\usage{
    \method{span}{TSestModel}(obj, compiled=.DSECOMPILED, ...)
}
\arguments{
    \item{obj}{A TSestModel }
}
\value{
The singular values of the matrix of tangent vectors are returned.
}
\sideffects{
If show.details is T then intermediate calculations are printed.
}
\details{
}
\seealso{
\code{\link{span}},
\code{\link{richardson.grad}},
\code{\link{genD}}
}

%\keyword{DSE}
\keyword{ts}
