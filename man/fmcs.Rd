\name{fmcs}
\alias{fmcs}
\title{
Flexible MCS (FMCS) Finder
}
\description{
R function to call the C++ implementation of the flexible common substructure (FMCS) algorithm. The FMCS algorithm provides an improved maximum common substructure (MCS) search method that allows atom and/or bond mismatches in the substructures shared among two small molecules. The resulting flexible MCSs (FMCSs) are often larger than strict MCSs, resulting in the identification of more common features in their source structures, as well as a higher sensitivity in detecting weak similarities among compounds.
}
\usage{
fmcs(sdf1, sdf2, al = 0, au = 0, bl = 0, bu = 0, matching.mode =
	  "static", fast = FALSE, timeout=0)
}
\arguments{
  \item{sdf1}{
Input query \code{SDF} object or \code{SDFset} object with a single molecule.
}
  \item{sdf2}{
Input target \code{SDF} object \code{SDFset} object with a single molecul.
}
  \item{al}{
Lower bound for the number of atom mismatches.
}
  \item{au}{
Upper bound for the number of atom mismatches.
}
  \item{bl}{
Lower bound for the number of bond mismatches.
}
  \item{bu}{
Upper bound for the number of bond mismatches.
}
  \item{matching.mode}{
Three modes for bond matching are supported: "static", "aromatic", and "ring".
}
  \item{fast}{
If \code{fast} is set to \code{TURE}, then the fast computing mode will be turned on. In this case, the algorithm will only return the size information about the source structures and their MCSs, while omitting all structural information.
}
	\item{timeout}{
		The maximum amount of time to spend searching, in milliseconds.
		A value of 0 indicates no timeout and is the default.
	}
}
\details{
...
}
\value{
Returns object of class \code{MCS}
}
\references{
Publication in preparation.
}
\author{
Yan Wang, Thomas Girke
}
\seealso{
\code{plotMCS}, \code{fmcsBatch}, \code{?"MCS-class"}
}
\examples{
library(fmcsR)
data(sdfsample)
sdfset <- sdfsample
mcs1 <- fmcs(sdfset[[1]], sdfset[[2]])
mcsfast <- fmcs(sdfset[[1]], sdfset[[2]], fast=TRUE)
mcs2 <- fmcs(sdfset[[1]], sdfset[[2]], au=3)
mcs3 <- fmcs(sdfset[[1]], sdfset[[2]], bu=3)
mcs4 <- fmcs(sdfset[[1]], sdfset[[2]], au=1, bu=1)
mcs5 <- fmcs(sdfset[[1]], sdfset[[2]], matching.mode="aromatic")
mcs6 <- fmcs(sdfset[[1]], sdfset[[2]], au=2, bu=1, matching.mode="aromatic")

## Plot MCS objects
plotMCS(mcs6)

## Methods to return components of MCS objects
stats(mcs6)
mcs6[["stats"]]
mcs1(mcs6)
mcs6[["mcs1"]]
mcs2(mcs6)
mcs6[["mcs2"]]

## Constructor method from list
mylist <- list(stats=stats(mcs6), mcs1=mcs1(mcs6), mcs2=mcs2(mcs6))
mymcs <- as(mylist, "MCS")
}
\keyword{ utilities }
