\name{fprint}
\alias{fprint}
\alias{fprint.tagged}
\title{Formatted Printing of Time Series}
\description{Generate a formatted character matrix of time series data.}
\usage{
    fprint(x, super.title=NULL, sub.title=NULL,
          digits=options()$digits, space=" ", file=NULL, append=FALSE)
    \method{fprint}{tagged}(x, super.title=NULL, sub.title=NULL, 
          digits=options()$digits, space=" ", file=NULL, append=FALSE) 
}
\arguments{
    \item{x}{a matrix of time series data.}
    \item{super.title}{a string to be used as top heading.}
    \item{sub.title}{a string to be used as second level of headings.}
    \item{digits}{number of digits to print.}
    \item{space}{a string to be used between columns of data.}
    \item{file}{name of a file to generate.}
    \item{append}{If TRUE output is appended to the file, otherwise 
      the file is overwritten.}
}
\value{A character matrix.}
\details{
    This is a generic method for print a table of time series data with
    formatting control. Currently the only specific method is fprint.tagged. 
}
\section{Side Effects}{If \code{file} is specified then a file will be created.}
\seealso{
    \code{\link{print}}
    \code{\link[tframe]{tfprint}}
}
\concept{DSE}
\keyword{ts}
\keyword{programming}
\keyword{utilities}

