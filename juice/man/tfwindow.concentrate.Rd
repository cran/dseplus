\name{tfwindow.concentrate}
\alias{tfwindow.concentrate}

\title{tfwindow Specific Methods}
\description{See the generic function description.}
\usage{
    \method{tfwindow}{concentrate}(x, tf=NULL,
              start=tfstart(tf), end=tfend(tf), warn=TRUE)
}
\arguments{
    \item{x}{a concentrate object to trucate.}
    \item{start}{A start date of a format compatible with the time series}
    \item{end}{An end date of a format compatible with the time series}
    \item{tf}{A tframe or tframed object}
    \item{warn}{A logical indicating if warning should be produced}
}

\seealso{ \code{\link[tframe]{tfwindow}} }

%\keyword{DSE}
\keyword{ts}

