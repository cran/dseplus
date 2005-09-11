\name{getpadi.default}
\alias{getpadi.default}

\title{Get Data from TSPADI Database Interface}
\description{Get data from a TSPADI database interface.}
\usage{
    \method{getpadi}{default}(series, server=Sys.info()[["nodename"]], dbname="",
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
 The function getpadi retrieves data from a (time series) server.
The padi code is also available (including the PADI server side)
at \url{www.bank-banque-canada.ca/pgilbert}.

(This documentation could use some cleaning up, but first I am really hoping to
 convert the whole interface to a more modern, non-RPC based, mechanism.)
  
start. and end dates or start dates and nobs can be supplied.  If all
are set to zero (the default) then all data is retrieved, provided
max.obs is large enough.  If more than one series is to be returned
then series should be a vector of strings. In this case dbname must be
a vector of corresponding length or a single string which is applied to
all series. All series must have the same frequency.  If the number of
observations in any single series is larger than max.obs then an error
will occur and max.obs should be set larger. For many data frequencies
startd and endd can be omitted.


 If start.server is FALSE then there will be no attempt to start a server and the 
      function will stop if a server is not running.
      If start.server is TRUE (the default) then if server==Sys.info()[["nodename"]] (the
      default) and a server is not already running, there will be an attempt
      to start a server using the argument server.process. The default is
      determined by PADIserverProcess()

 The argument server.process is only used if it is necessary to start a server.

 If specified, server would typically be a single string, though there is
    some attempt to handle vectors of strings (indicating different servers
    for each series).

 The string dbname is passed to the server, but it may or may not be used, 
    depending on the server implementation.
    If the server is being started then it will run in the Unix pwd and local
    path names should work, BUT in general there is no guarantee that the 
    server is running in the pwd and complete path names may be required.

 Start and end dates or start dates and nobs can be supplied.
 If all are set to zero (the default) then all data is retrieved, provided 
    max.obs is large enough.

 If the starting date and ending date and number of observations is set to
 zero, it will return the whole series.  If the starting and ending dates
 are zero but number of observations is set, the LAST numobs observations
 are returned.  If the starting date is sent, ending date is zero and numobs
 is non-zero, the FIRST numobs starting from start date are returned.

 The size of the data array passed to C is the larger of nobs and max.obs.
 A ts matrix is returned. If more than one series is to be returned then
  series should be a vector of strings. In this case server must
  be a vector of corresponding length or a single string which is applied to all
  series. All series must have the same frequency and the time window
  is the intersection of the window for each series (i.e. the latest start and
  earliest end date).
 startd and endd provide for a tag (day).

 If the starting date and ending date and number of observations is set to
 zero, it will return the whole series.  If the starting and ending dates
 are zero but number of observations is set, the LAST numobs observations
 are returned.  If the starting date is sent, ending date is zero and numobs
 is non-zero, the FIRST numobs starting from start date are returned.

 If transformations is not null it should be a vector of strings, one for 
 each series (with "" for any series which is not to be transformed), or a
 single string which is applied to each series. 
 The transformations are applied by eval(call(transformations[i], .).

 If use.tframe=FALSE then ts() is used to construct the time series, otherwise
    the tframe approach is used. (See the tframe package.)
    Certain transformations available with DSE require the tframe structure
    and an error may result if these transformations are attempted with
    use.tframe=FALSE.

 If stop.on.error is TRUE then stop is executed if a data retrieval error
 occurs. If stop.on.error is FALSE then the error message is returned, which means
 the calling function will need to handle the error.

 If warn is TRUE then warnings are issued for certain data conversions (weekly 
 data) which may not work in all case. If warn is these warnings are not issued.
 timeout is the period to wait (in seconds) before generating an error.
}
\seealso{
   \code{\link{checkPADIserver}}
   \code{\link{putpadi.default}}
   \code{\link[dsepadi]{getpadi}}
   \code{\link[dsepadi]{putpadi}}
}
\examples{
if(require("padi") && checkPADIserver("ets")){
  cpi <-getpadi("P100000", server= "ets")# June 1992=100
  cpi <-getpadi("P100000", starty=1988, startm=1,endy=1990, endm=12, server= "ets")
  cpi <-getpadi( "P100000",starty=1988, startm=1,nobs=8, server= "ets")
  data <- getpadi( c("V122491","V37124","P100000"), server= "ets" )

  cpi <-getpadi("P100000", server= "ets", use.tframe=TRUE)
  data <- getpadi( c("V122491","V37124","P100000"), server= "ets", use.tframe=TRUE)
 }
}
%\keyword{tframe}
\keyword{ts}


