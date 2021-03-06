\name{setTSPADIdata}
\alias{setTSPADIdata}
\title{Prompt for Series Names}
\description{Prompt for data series names and locations and create a data object.}
\usage{
    setTSPADIdata()
}
\value{An object of class 'TSdata' or of class 'TSPADIdata' 'TSdata'}
\details{
The function prompts for the names and database locations of series
to be used for input (exogenous) variables and output (endogenous)
variables and creates a data object. The data can be retrieved, in 
which case it is fixed as currently available (and becomes an object
of class 'TSdata'), or it may be left dynamic, in which case it is
retrieved whenever it is referenced (and is of class 'TSPADIdata' 'TSdata').
}
\seealso{
\code{\link[tframe]{freeze}}
\code{\link[dse1]{TSdata}}
}
\examples{
\dontrun{data <- setTSPADIdata()}
}
%\keyword{DSE}
\keyword{ts}

