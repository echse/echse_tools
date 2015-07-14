\name{geogrid.readAscii}
\alias{geogrid.readAscii}
\title{Read a geogrid object from a file}
\usage{
  geogrid.readAscii(file)
}
\arguments{
  \item{file}{Name/path of the ASCII grid file.}
}
\value{
  An object of class geogrid.
}
\description{
  The expected file format is ESRI's ASCII grid format,
  i.e. there must be six header lines with keywords
  'ncols', 'nrows', 'xllcorner', 'yllcorner', 'cellsize',
  and 'nodata_value'. The matrix of values starts at line
  seven.
}
\examples{
# Create a sample file
g= geogrid(matrix(runif(6*8), ncol=8, nrow=6, byrow=TRUE),
  xllcorner=0, yllcorner=0, cellsize=1, nodata=-9999)
f= tempfile()
g= geogrid.writeAscii(grid=g, file=f)
# Read and display
g= geogrid.readAscii(f)
geogrid.plot(g,"Random data")
file.remove(f)
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}

