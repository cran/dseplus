\name{getfm}
\alias{getfm}
\alias{putfm}

\title{Function is for backwards compatability at the BOC}
\description{Function is for backwards compatability at the BOC}
\usage{
    getfm(dbname,series, starty=0,startp=0, endy=0,endp=0, 
       nobs=0,max.obs=2000, transformations=NULL, pad=FALSE)
    putfm(data, dbname, seriesNames)
}
\arguments{
    \item{dbname} {a string indicating Fame database name.}

    \item{series}{see \code{getpadi}.}
    \item{starty}{see \code{getpadi}.}
    \item{startp}{start period.}
    \item{endy}{see \code{getpadi}.}
    \item{endp}{end period.}
    \item{nobs}{see \code{getpadi}.}
    \item{max.obs}{see \code{getpadi}.}
    \item{transformations}{see \code{getpadi}.}
    \item{pad}{see \code{getpadi}.}
    \item{data}{a object with data to put on the server.}
    \item{seriesNames}{vector of strings to use for series identifiers.}
}
\value{{see \code{getpadi}.}}
\details{
   Depreciated.
}
\seealso{
\code{\link{getpadi}}

}
%\examples{
%}
%\keyword{DSE}
\keyword{ts}

