\name{plot2by2}
\alias{plot2by2}
\alias{plot2by2.default}
\alias{plot2by2.TSdata}

\title{plot2by2}
\description{
    plot data series one vs another, two at a time (that is, data[,i] vs 
    data[,j] for all i,j (not on time axis).
}
\usage{
    plot2by2(data, ...)
    \method{plot2by2}{default}(data, pch=".", ...)
    \method{plot2by2}{TSdata}(data, ...)
}
\arguments{
    \item{data}{a matrix of time series or a TSdata object.}
    \item{pch}{character to be used for plotting.}
    \item{...}{arguments passed to tfplot.}
}
\value{None.}
\sideaffects{A plot is produced.}

\concept{DSE}
\keyword{ts}

