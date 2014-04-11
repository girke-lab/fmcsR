\name{fmcsR-package}
\alias{fmcsR-package}
\alias{fmcsR}
\docType{package}
\title{
%What the package does (short line)
%~~ package title ~~
A FMCS solver package.
}
\description{
%More about what it does (maybe more than one line)
%~~ A concise (1-5 lines) description of the package ~~
The package consists of two main functions, \code{fmcs} which computes the flexible MCS between two \code{SDF} objects. And \code{fmcsBatch} runs the FMCS algorithm on a \code{SDFset}.
}
\details{
\tabular{ll}{
Package: \tab fmcsR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-02-01\cr
}
%~~ An overview of how to use the package, including the most ~~
%~~ important functions ~~
}
\author{
Yan Wang

Maintainer: Yan Wang <wangya@cs.ucr.edu>

}

%~~ Optionally other standard keywords, one per line, from file ~~
%~~ KEYWORDS in the R documentation directory ~~
\keyword{ package }

\examples{

library(fmcsR)
data(sdfsample)
sdfset <- sdfsample
result1 <- fmcs(sdfset[[1]], sdfset[[2]])
result2 <- fmcs(sdfset[[1]], sdfset[[2]], au=3)
result3 <- fmcs(sdfset[[1]], sdfset[[2]], bu=3)
result4 <- fmcs(sdfset[[1]], sdfset[[2]], au=1, bu=1)
result5 <- fmcs(sdfset[[1]], sdfset[[2]], matching.mode="aromatic")
result6 <- fmcs(sdfset[[1]], sdfset[[2]], au=2, bu=1, matching.mode="aromatic")

fmcsBatch(sdfset[[1]], sdfset[1:3])
fmcsBatch(sdfset[[1]], sdfset[1:3], au=2)
fmcsBatch(sdfset[[1]], sdfset[1:3], bu=1)
fmcsBatch(sdfset[[1]], sdfset[1:3], matching.mode="aromatic", au=1, bu=1)

}
