\name{tags}
\alias{tags}
\alias{tags<-}
\alias{tagged}
\alias{tagged.default}
\alias{tagged.TSdata}
\alias{is.tagged}

\title{Tagged Matrices}
\description{Matrices with an decription for each point. }
\usage{
    tags(x)
    tags(x) <- value
    tagged(x, ...)
    \method{tagged}{default}(x, tags)
    \method{tagged}{TSdata}(x, input.tags, output.tags)
    is.tagged(obj)
}
\arguments{
    \item{x}{a matrix or TSdata object.}
    \item{value}{if x is a matrix then \code{value} should be a matrix of 
       strings of the same dimension or a scalar string which is expanded
        to the dimension of x.}
    \item{tags}{if x is a matrix then \code{tags} should be a matrix of 
       strings of the same dimension or a scalar string which is expanded
        to the dimension of x.}
    \item{input.tags}{
        If x is TSdata then input.value should be a matrix of strings of 
        the same dimension as input.data(x)or a single string which is 
        expanded.}
    \item{output.tags}{
        If x is TSdata then output.value should be a matrix of strings of 
        the same dimension as input.data(x)or a single string which is 
        expanded.}
    \item{...}{arguments to be passed to other methods.}
    \item{obj}{any object.}
}
\value{
    A matrix of class "tagged" or a TSdata object with specific 
    class "tagged"}
\details{
    The matrix of descriptive information is assigned as an attribute of the
    matrix. Most testing of the methods for this class has been with a 
    single character tag which can be used as a flag, for example, to 
    indicate the different sources for the data points.
}
\seealso{
    \code{\link{simple.monitoring}}
    \code{\link{TSdata}}
}
\examples{
    x <- tagged(matrix(rnorm(100), 50,2), "r")
    is.tagged(x)
}

%\keyword{DSE}
\keyword{ts}

