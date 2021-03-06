\name{00.juice.Intro}
\alias{00.juice.Intro}
\alias{juice.Intro}

\docType{package}

\title{Juice}
\description{Methods for Concentrating and Reconstituting Data}
\details{
    The attempted approach is that a model of type TSmodelconcentrate  
    should work like KF and ARMA model, so that plot, residual, 
    etc., produce results for the full (reconstitued) data set. 
    Special methods (eg. concentrated.tfplot,
    concentrated.checkResiduals) do the equivalent thing 
    using the reduced dimensional data as the TSdata and 
    dropping the fact that the model is of 
    class TSmodelconcentrate. 

     The outputData for a concentrated object is the original 
    data, but for a reconstituted object it is the reconstituted data. 
    So, for example, \code{tfplot(ConcentratedDataObject)} will plot 
    the original data (as well as the reconstructed data) whereas 
    \code{tfplot(reconstitute(ConcentratedDataObject))}
    plots only the reconstituted data. (To plot the concentrated 
    series use \code{concentrated.tfplot(ConcentratedDataObject)}. 

    Typically one should not work with a reconstituted object unless it is 
    explicitely needed, as the original data is "hidden".

    The concentrated data is extracted from both 
    \code{TSdataconcentrate} and \code{TSdatareconstitute} using 
    the function \code{concentrated.outputData}.
}
\concept{DSE}
\keyword{ts}

