\name{tfwindow.tagged}
\alias{tfwindow.tagged}
\alias{tbind.tagged}
\alias{trim.na.tagged}
\alias{splice.tagged}
\alias{select.series.tagged}

\title{Specific Methods for tframed Data}
\description{See the generic function description.}
\usage{
    \method{tfwindow}{tagged}(x, start.=NULL, end.=NULL, tf=NULL, warn=T)
    \method{tbind}{tagged}(mat1, mat2)
    \method{trim.na}{tagged}(mat, start.=T, end.=T)
    \method{splice}{tagged}(mat1, mat2, tag1=tags(mat1), tag2=tags(mat2))
    \method{select.series}{tagged}(x, series=seq(ncol(x)))
}
\arguments{
    \item{x}{see the generic function.}
    \item{start.}{see the generic function.}
    \item{end.}{see the generic function.}
    \item{tf}{see the generic function.}
    \item{warn}{see the generic function.}
    \item{mat}{see the generic function.}
    \item{mat1}{see the generic function.}
    \item{mat2}{see the generic function.}
    \item{tag1}{tags for the first matrix. See tags.}
    \item{tag2}{tags for the first matrix. See tags.}
    \item{series}{see the generic function.}
}
\seealso{
    \code{\link{tfwindow}},
    \code{\link{tbind}},
    \code{\link{trim.na}}
    \code{\link{splice}}
    \code{\link{select.series}}
}

%\keyword{DSE}
\keyword{ts}

