\name{tframed.concentrate}
\alias{tframed.concentrate}
\alias{tframe<-.concentrate}

\title{Construct a Tframed Object}
\description{Create a trfamed object or set the tframe of an object.}

\usage{
    \method{tframed}{concentrate}(x, tf=NULL, names = NULL)
    \method{tframe}{concentrate}(x) <- value
}

\arguments{
    \item{x}{a (tframed) concentrate object or a concentrate object to be tframed.}
    \item{tf}{a tframe attribute to be applied to x.}
    \item{names}{optional (new) series names to be applied to x.}
}
\value{A tframed object.}
\details{See the generic.
}
\seealso{
    \code{\link[tframe]{tframe}}
}

\keyword{programming}
\keyword{utilities}
\keyword{internal}
