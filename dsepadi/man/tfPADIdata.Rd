\name{tfPADIdata}
\alias{tfPADIdata}
\alias{is.tfPADIdata}
\alias{[.tfPADIdata}

\title{make tfPADIdata object}
\description{Make a tfPADIdata object from a vector of identifiers}
\usage{
    tfPADIdata(series,  server = "", db= "", transforms= "",  
           start=NA, end=NA, frequency=NA, names=NULL, 
           pad=FALSE, pad.start=pad, pad.end=pad,
           use.tframe=TRUE,
           start.server=NULL, server.process=NULL, cleanup.script=NULL,
           stop.on.error=TRUE, warn=TRUE)
    is.tfPADIdata(obj)
}
\arguments{
    \item{series}{vectors of strings indicating series identifiers.}
    \item{server}{string or vector of strings indicating the server where
       the series will be found. Scalar values are expanded to a 
       vector of appropriate length.}
    \item{db}{string or vector of strings indicating additional information for 
       the server. Scalar values are expanded to a vector of appropriate length.}
    \item{transforms}{see details.}
    \item{start}{a two element vector of starting year and period which 
       is used when requesting data.}
    \item{end}{a two element vector of end year and end period which 
       is used when requesting data.}
    \item{frequency}{a scalar indicating frequency. This is not 
       used in requesting data, but if supplied it is used as a check of 
       returned data. A warning is issued if returned data does not have this 
       frequency.}
    \item{names}{names to replace series identifiers.}
    \item{pad}{TRUE or FALSE indicating if data should be padded with NAs on both ends
       to the length of the longest available series. If padding is not done 
       then series are truncated to the intersection of available data so there 
       are no NAs in the result. pad.start and pad.end are the arguments which are 
       actually used. pad is just a convenient way to specify both.}
    \item{pad.start}{logical indicating if NAs should be padded to beginning of data.}
    \item{pad.end}{logical indicating if NAs should be padded to end of data.}
    \item{use.tframe}{logical indicating if the tframe library methods should be
       used for the time dimension.}
    \item{start.server}{logical indicating if a (local) server should be
       started.}
    \item{server.process}{string indicating the name to be used to start a
       server process.}
    \item{cleanup.script}{string indicating the name to be used to shut down a
       server process.}
    \item{stop.on.error}{logical indicating if the function should stop if any
       series produces an error, or continue with other series.}
    \item{warn}{logical indicating if warning messages should be supressed.}
    \item{obj}{Any object.}
}

\value{
A object of class "tfPADIdata" which can be used to retrieve a matrix 
time series object.}
\details{
This is the way to produce an object which can be used for database access. 
The transforms a vector of strings, one for each series, giving any 
transformations which should be applied when the data is retrieved (with eval(call()), e.g. "log",  or "" if no transform is to be applied).
}
\seealso{
\code{\link{freeze}}
}
\examples{
    z <- tfPADIdata( c("seriesA", "seriesB"), server="myserver")
    is.tfPADIdata(z)
}
%\keyword{tframe}
\keyword{ts}

