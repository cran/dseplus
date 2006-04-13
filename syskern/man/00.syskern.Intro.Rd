\name{00.syskern.Intro}
\alias{00.syskern.Intro}
\alias{syskern.Intro}

\docType{package}

\title{Syskern}

\section{Usage}{
   library("syskern") 
}

\description{
  A small set of programs to facilitate programming around S/R differences}
\section{Introduction}{
The purpose of this package is to provide a 
common set of utilities for some small points where S and R differ, 
so that other code works in both S and R, and to some
extent on different operating systems. Since R now supports different operating
systems fairly well the package has been revised to use R function 
names and syntax. A skeleton of similar programs is provided for use in S
(tested only with Splus 3.3 in Unix). 

The files in the S/ subdirectory of the package are for use with S (not R) and
provide the functions \code{is.R}, \code{is.S}, \code{Sys.getenv}, \code{Platform}, 
\code{Sys.sleep}, \code{Sys.info}, \code{Sys.mail}, \code{system}, 
\code{file.copy}, \code{getwd}, \code{getOption}. The last of these uses the
global variable \code{.Syskern.options} which is used only to indicate
the mailer option. In many cases these functions just call S equivalents, 
possibly with arguments different from the defaults. The function \code{require} 
is defined but does nothing in S.

There are also a few function and variables added which have different
definitions in S and R. A number of \code{is.Sxx} functions are defined to 
identify versions of S. The variable \code{.SPAWN} is also set to be used in
Splus to identify if "For" loops should be used. (It is sometimes better not 
to use these even in Splus. \code{.SPAWN} should be FALSE in R.)

%The default action of \code{unlink} is to remove subdirectories 
%in the Unix version of Splus 3.3 and
%this cannot be reconciled with R which requires the non-default action to 
%remove subdirectories. Thus it has been necessary to define the
%function \code{syskern.rm} which recursively removes subdirectories.

The functions code{file.date.info}, \code{dateParsed} used in the dse bundle
are also defined but are being depreciated and should be avoided.
}
  
\keyword{programming}
\keyword{interface}
\keyword{utilities}
\keyword{distribution}

\author{Paul Gilbert}
