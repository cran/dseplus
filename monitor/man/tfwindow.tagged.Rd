\name{tfwindow.tagged}
\alias{tfwindow.tagged}
\alias{tbind.tagged}
\alias{splice.tagged}
\alias{selectSeries.tagged}
\alias{tframe<-.tagged}

\title{Specific Methods for tframed Data}
\description{See the generic function description.}
\usage{
    \method{tfwindow}{tagged}(x, tf=NULL, start=tfstart(tf), end=tfend(tf), warn=TRUE)
    \method{tbind}{tagged}(x, mat2, ..., pad.start=TRUE, pad.end=TRUE, warn=TRUE)
    \method{splice}{tagged}(mat1, mat2, tag1=tags(mat1), tag2=tags(mat2), ...)
    \method{selectSeries}{tagged}(x, series=seq(ncol(x)))
    \method{tframe}{tagged}(x) <- value
}
\arguments{
    \item{x}{see the generic function.}
    \item{start}{see the generic function.}
    \item{end}{see the generic function.}
    \item{tf}{see the generic function.}
    \item{pad.start}{see the generic function.}
    \item{pad.end}{see the generic function.}
    \item{warn}{see the generic function.}
    \item{mat1}{see the generic function.}
    \item{mat2}{see the generic function.}
    \item{tag1}{tags for the first matrix. See \code{tags}.}
    \item{tag2}{tags for the second matrix. See \code{tags}.}
    \item{series}{see the generic function.}
    \item{value}{see the generic function.}
    \item{...}{(further arguments, currently disregarded)}
}
\seealso{
    \code{\link[tframe]{tfwindow}},
    \code{\link[tframe]{tbind}},
    \code{\link[tframe]{trimNA}}
    \code{\link[tframe]{splice}}
    \code{\link[tframe]{selectSeries}}
}

\concept{DSE}
\keyword{ts}

