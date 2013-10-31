\name{MCS-class}
\Rdversion{1.1}
\docType{class}
\alias{MCS-class}
\alias{[[,MCS-method}
\alias{coerce,list,MCS-method}
\alias{mcs1,MCS-method}
\alias{mcs1}
\alias{mcs2,MCS-method}
\alias{mcs2}
\alias{stats,MCS-method}
\alias{stats}
\alias{show,MCS-method}

\title{Class \code{"MCS"}}
\description{
List-like container for storing results from \code{fmcs} function. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MCS", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{stats}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{mcs1}:}{Object of class \code{"SDFset"} ~~ }
    \item{\code{mcs2}:}{Object of class \code{"SDFset"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{[[}{\code{signature(x = "MCS")}: ... }
    \item{coerce}{\code{signature(from = "list", to = "MCS")}: ... }
    \item{mcs1}{\code{signature(x = "MCS")}: ... }
    \item{mcs2}{\code{signature(x = "MCS")}: ... }
    \item{stats}{\code{signature(x = "MCS")}: ... }
	 }
}
\references{
...
}
\author{
Yan Wang
}
\note{
...
}

\seealso{
Related classes: SDF, SDFstr
}
\examples{
## Create MCS instance
showClass("MCS")
data(sdfsample)
sdfset <- sdfsample
mcs <- fmcs(sdfset[[1]], sdfset[[2]], au=2, bu=2)

## Methods to return components of MCS 
stats(mcs)
mcs[["stats"]]
mcs1(mcs)
mcs[["mcs1"]]
mcs2(mcs)
mcs[["mcs2"]]

## Constructor method from list
mylist <- list(stats=stats(mcs), mcs1=mcs1(mcs), mcs2=mcs2(mcs))
mymcs <- as(mylist, "MCS")
}
\keyword{classes}
