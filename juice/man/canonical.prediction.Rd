\name{canonical.prediction}
\alias{canonical.prediction}
\alias{is.TScanonicalPrediction}

\title{Canonical Prediction}
\description{
    Use canonical correlation with input data as the independent variables 
    used to predict output data.}
\usage{
    canonical.prediction(d, conc=concentrator(d),
       q=min(concentrated.nseriesInput(d),
             concentrated.nseriesOutput(d)))
    is.TScanonicalPrediction(x)
}
\arguments{
    \item{d}{a TSdataconcentrate object as returned by concentrate.}
    \item{conc}{a concentrator.}
    \item{q}{integer indicating the number of canonical variates to keep.}
    \item{x}{any object.}
}
\value{A TScanonicalPrediction object.}
\details{
    Data d as returned by concentrate. Alternately, a different
    conc (proj) can be used.
    Use q canonical variates from  input data as predictors of 
        q canonical variates from output data and then use these to 
       reconstruct output data. (ref T.W. Anderson p491)
    q cannot exceed 
     min(concentrated.nseriesInput(d), concentrated.nseriesOutput(d)) 
}
\seealso{
    \code{\link{concentrate}}
    \code{\link{concentrator}}
}

\examples{
    if(is.R()) data("eg1.DSE.data.diff", package="dse1")
    z <- canonical.prediction(concentrate(eg1.DSE.data.diff))
    is.TScanonicalPrediction(z)
}
\concept{DSE}
\keyword{ts}


