\name{input.data.TSPADIdata}
\alias{input.data.TSPADIdata}
\alias{periodsInput.TSPADIdata}
\alias{output.data.TSPADIdata}
\alias{periodsOutput.TSPADIdata}

\title{TS Input and Output Specific Methods}
\description{See the generic function description.}
\usage{
    \method{input.data}{TSPADIdata}(x, series=seq(length=nseriesInput(x)))
    \method{output.data}{TSPADIdata}(x,series=seq(length=nseriesOutput(x)))
    \method{periodsInput}{TSPADIdata}(x)
    \method{periodsOutput}{TSPADIdata}(x)
}
\arguments{
    \item{x}{a TSPADIdata object.}
    \item{series}{series to select. Passed to select.series.}
}
%\keyword{DSE}
\keyword{ts}


