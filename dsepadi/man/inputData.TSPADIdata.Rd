\name{inputData.TSPADIdata}
\alias{inputData.TSPADIdata}
\alias{periodsInput.TSPADIdata}
\alias{outputData.TSPADIdata}
\alias{periodsOutput.TSPADIdata}

\title{TS Input and Output Specific Methods}
\description{See the generic function description.}
\usage{
    \method{inputData}{TSPADIdata}(x, series=seq(length=nseriesInput(x)))
    \method{outputData}{TSPADIdata}(x,series=seq(length=nseriesOutput(x)))
    \method{periodsInput}{TSPADIdata}(x)
    \method{periodsOutput}{TSPADIdata}(x)
}
\arguments{
    \item{x}{a TSPADIdata object.}
    \item{series}{series to select. Passed to selectSeries.}
}
%\keyword{DSE}
\keyword{ts}


