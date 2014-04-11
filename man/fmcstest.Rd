\name{fmcstest}
\alias{fmcstest}
\docType{data}
\title{
SD file stored in \code{SDFset} object
}
\description{
Sample compound structures stored in SDF format.
}
\usage{data(fmcstest)}
\format{
  Object of class \code{SDFset}
}
\details{
Object stores X molecules from a sample SD file.
}
\source{
ftp://ftp.ncbi.nih.gov/pubchem/
}
\references{
SDF format definition: http://www.symyx.com/downloads/public/ctfile/ctfile.jsp
}
\examples{
data(fmcstest)
sdfset <- fmcstest
view(sdfset)
}
\keyword{datasets}
