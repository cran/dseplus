\name{simple.monitoring}
\alias{simple.monitoring}
\alias{check.for.value.changes}
\alias{check.for.file.date.changes}
\alias{watch.data}
 
\title{Simple Monitoring}
\description{Automatic monitoring with e-mail of results.}
\usage{
    simple.monitoring(model, data.names, 
   previous.data=NULL,
   mail.list=NULL,
   error.mail.list=Sys.info()[["user"]],
   message.title="Simple Monitoring",
   message.subject="Simple Monitoring",
   message.footnote=NULL,
   show.start= c(0,-3),
   show.end  = c(0,12),    
   report.variables= seriesNames(data.names),
   data.sub.heading=NULL,
   data.tag=" ",
   forecast.tag="f",
   run.again=FALSE,
   save.as=NULL)
}
\arguments{
    \item{model}{a TSmodel.}
    \item{data.names}{a TSdata (names) object.}
    \item{previous.data}{a TSdata object similar to that returned by
       freeze(data.names). Used to check if the data has been updated.}
    \item{mail.list}{string containing user ids for mail}
    \item{error.mail.list}{string contain user ids for mail}
    \item{message.title}{string or vectors of strings placed at the top 
       of the message.}
    \item{message.subject}{string used as the message subject.}
    \item{message.footnote}{string or vectors of strings placed at the bottom 
       of the message.}
    \item{show.end}{integer indicating the number of periods after the 
       end of data which should be displayed.(i.e. number of forecast periods)}
    \item{show.start}{negative integer indicating the number of periods 
       before the end of data (i.e. history) which should be displayed. It 
       is added to the end (so it should be negative).}
    \item{report.variables}{indicates output variables which are reported.
       It should be a vector of strings corresponding to a subset of 
       names returned by seriesNames.}
    \item{data.sub.heading}{string or vectors of strings placed at the top 
       of the data.}
    \item{data.tag}{matrix of strings placed beside data points in 
       the report.}
    \item{forecast.tag}{matrix of strings placed beside data points in 
       the report.}
    \item{run.again}{logical indicating that the monitoring should be run 
       regardless of data updates.}
    \item{save.as}{optional string giving file name in which to save details of
       the data and model (useful for debugging).}
 }
\value{
    Invisibly return latest data for next comparison. This function is
    run mainly for its side effects.}
\details{
    \code{mail.list} and \code{error.mail.list} should be single strings (not 
    vectors) but the string can contain multiple user ids for mail.
    If mail.list is NULL (default) then mail is not sent (useful for testing). 
    If error.mail.list  is  NULL   then mail is not sent (useful for testing).
    The default for error.mail.list is the result of Sys.info()[["user"]].
    This version does not allow for 
        -combining forecasts (ie. monitoring.data or overriding data)
        -input (policy) projections
    See \code{combination.monitoring} for these features.
}
\sideffects{Mail is sent.}
\seealso{
    \code{\link{combination.monitoring}}
    \code{\link[syskern]{Sys.mail}}
}

%\keyword{DSE}
\keyword{ts}

