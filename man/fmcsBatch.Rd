\name{fmcsBatch}
\alias{fmcsBatch}
\title{
FMCS Search Function
}
\description{
Compound search function that runs the FMCS algorithm for a query compound against a set of molecules stored
in an \code{SDFset} container.
}
\usage{
fmcsBatch(querySdf, sdfset, al = 0, au = 0, bl = 0, bu = 0,
matching.mode = "static",timeout=0)
}
\arguments{
  \item{querySdf}{
Input query \code{SDF} object or \code{SDFset} object of length one.
}
  \item{sdfset}{
Input target \code{SDFset} object.
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
Three matching mode are supported, "static", "aromatic", and "ring".
}
	\item{timeout}{
		The maximum amount of time to spend on each pair of
		comparisons, in milliseconds. The default value of 0 indicates no timeout.
	}
}
\details{
This function runs the FMCS algorithm in fast computing mode. Thus, it will only return the similarity scores and size information about the source structures and their MCSs, while omitting all structural information.
}
\value{
Returns a \code{matrix} with compound IDs as row names and the following columns:
Query_Size, Target_Size, MCS_Size, Tanimoto_Coefficient and
Overlap_Coefficient. For details see vignette of this package.
}
\author{
Yan Wang, Thomas Girke
}
\seealso{
\code{plotMCS}, \code{fmcs}, \code{?"MCS-class"}
}
\examples{
library(fmcsR)
data(sdfsample)
sdfset <- sdfsample
fmcsBatch(sdfset[[1]], sdfset[1:3])
fmcsBatch(sdfset[[1]], sdfset[1:3], au=2)
fmcsBatch(sdfset[[1]], sdfset[1:3], bu=1)
fmcsBatch(sdfset[[1]], sdfset[1:3], matching.mode="aromatic", au=1, bu=1)
}
\keyword{utility}
