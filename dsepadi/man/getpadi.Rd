\name{getpadi}
\alias{getpadi}

\title{Get Data from TSPADI Database Interface}
\description{Get data from a TSPADI database interface.}
\usage{
    getpadi(series, server=Sys.info()[["nodename"]], dbname="",
                start.server=TRUE, server.process=PADIserverProcess(),
                cleanup.script=PADIcleanupScript(),
                starty=0,startm=0,startd=1, endy=0,endm=0,endd=1, 
                nobs=0,max.obs=2000, transformations=NULL, pad=FALSE,
                user=Sys.info()[["user"]], passwd="",
                stop.on.error=TRUE, use.tframe=FALSE, warn=TRUE, timeout=60)

}
\arguments{
    \item{series}{
       A character string giving the name of the series. Alternately, 
       series can be a vector of character strings specifying multiple series.}
    \item{server}{
       A character string giving the network name of the server which is to be 
       requested to supply the series.  If series specifies multiple series and 
       they are not all on the same server then server should be a vector of 
       character strings with elements corresponding to the elements of series.}
    \item{dbname}{
       A character string giving additional information to the server about the 
       location of the series (eg. the name of a database).  If series specifies 
       multiple series and they are not all on the same database then dbname 
       should be a vector of character strings with elements corresponding to 
       the elements of series.}
    \item{starty}{An integer indicating the starting year.}
    \item{startm}{An integer indicating the starting period.}
    \item{startd}{An integer indicating the starting day.}
    \item{endy}{An integer indicating the ending year.}
    \item{endm}{An integer indicating the ending period.}
    \item{endd}{An integer indicating the ending day.}
    \item{nobs}{The number of observations.}
    \item{max.obs}{integer indicating the possible returned data size used
        to define the size of the buffer to prepare.}
    \item{transformations}{
       A character string giving transformations to be applied to the series 
       (e.g. "log"). If multiple series are being requested then transformations 
       can be a single string, in which case it is applied to all series, or a
       vector of character strings, one for each series. If no transformation is to
       be applied to some series then "" should be used.}
    \item{pad}{
       If FALSE (default) then all series are truncated to the interection of
       available time periods (i.e. the latest start date and earliest end date).
       If TRUE then series are padded with NA so the result
       starts at the earliest available observation and ends at the last available
       observation.}
    \item{start.server}{try to start a server if one is not running}
    \item{server.process}{command to execute in an attempt to start a server}
    \item{cleanup.script}{command to execute to terminate a server if one is started}
    \item{user}{user id for access to the database (if necessary)}
    \item{passwd}{password for access to the database (if necessary)}
    \item{stop.on.error}{
       If TRUE then stop is executed when an error occurs. Otherwise, the
       error message is returned and the calling program must deal with it.}
    \item{use.tframe}{
       If use.tframe=FALSE then ts() is used to construct the time series, otherwise
       the tframe utilities are used.
       Certain transformations available with DSE require the tframe stucture
       and an error may result if these transformations are attempted with
       use.tframe=FALSE.}
    \item{warn }{
       Print warning messages for some crude frequency conversions (weekly data).}
    \item{timeout}{an integer indicating the number of seconds to wait before
       concluding that the server is not available.}
}
\value{A time series matrix with a column for each series.}
\details{
   See getpadi.default}
\seealso{
%%   \code{\link{checkPADIserver}}
   \code{\link{putpadi}}
%%   \code{\link[padi]{getpadi.default}}
}

%\keyword{tframe}
\keyword{ts}


