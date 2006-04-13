\name{concentrated.tfplot}
\alias{concentrated.tfplot}
\title{Plot Concentrated Series}
\description{
    The concentrate data is plotted.
}
\usage{
    concentrated.tfplot(x, ...)
}
\arguments{
    \item{x}{A concentrated data object.}
    \item{...}{arguments to be passed to other tfplot.}
}
\value{
    Depends on the argument. For a simple concentrated data object 
    the result is a vector of strings.
}
\section{Side Effects}{A plot is generated}

\seealso{
    \code{\link[tframe]{tfplot}}
}
\examples{
    data("eg1.DSE.data", package="dse1")
    require("stats")
    z <- concentrate(eg1.DSE.data)
    concentrated.tfplot(z)
}
\concept{DSE}
\keyword{ts}

