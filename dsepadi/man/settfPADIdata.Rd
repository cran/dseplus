\name{settfPADIdata}
\alias{settfPADIdata}
\title{Prompt for series names}

\description{Prompt for data series names and locations and 
create a data object.}
\usage{settfPADIdata(preamble=TRUE)}
\arguments{
    \item{preamble}{a logical indicating if the description of how to enter data
    should be printed.}
}
\value{An object of class 'tfPADIdata' (or retrieved time series matrix).}
\details{
The function prompts for the names and database locations of series.
The data can be retrieved, in which case it is returned as currently 
available on the database (and becomes a time series matrix), or it may 
be left dynamic, in which case it is of class 'tfPADIdata' and the discription can be used anytime to retrieve data.
}
\seealso{
\code{\link[tframe]{freeze}}
\code{\link[dse1]{TSdata}}
\code{\link[dse1]{TSdata.object}}
}
\examples{
  data <- settfPADIdata()
 }
%\keyword{tframe}
\keyword{ts}

