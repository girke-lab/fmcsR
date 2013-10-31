\name{mcs2sdfset}
\alias{mcs2sdfset}
\title{
Return MCS object as SDFset
}
\description{
Helper function to run \code{atomsubset} from \code{ChemmineR} library on MCS
objects in order to obtain their results in \code{SDFset} format. }
\usage{
mcs2sdfset(x, ...)
}
\arguments{
  \item{x}{
Object of class \code{MCS}
}
  \item{\dots}{
Arguments to be passed to/from other methods.
}
}
\details{
Returns MCS data in form of a list containing two \code{SDFset} objects, one for the query and one for the target structure.
}
\value{
List with two \code{SDFset} objects.
}
\references{
...
}
\author{
Thomas Girke
}
\note{
...
}
\seealso{
\code{fmcs}
}

\examples{
library(fmcsR)
data(sdfsample)
sdfset <- sdfsample
mcs <- fmcs(sdfset[[1]], sdfset[[2]], au=2, bu=1, matching.mode="aromatic")
mcs2sdfset(x=mcs, type="new")
mcs2sdfset(x=mcs, type="old")[[1]][[1]] 
plot(mcs2sdfset(x=mcs, type="new")[[1]][1]) 
}
\keyword{ utilities }

