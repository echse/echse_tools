\name{northangle}
\alias{northangle}
\title{Point of the compass when looking from one point to another point}
\usage{
  northangle(x_from, y_from, x_to, y_to)
}
\arguments{
  \item{x_from}{x-coordinate of the base point}

  \item{y_from}{y-coordinate of the base point}

  \item{x_to}{x-coordinate of the point located from the
  base point}

  \item{y_to}{y-coordinate of the point located from the
  base point}
}
\value{
  An angle in range 0...360.
}
\description{
  When looking from a point 'A' to another point 'B', the
  function computes the point of the compass corresponding
  to 'B'. The result is the angle against North (0 ==
  North, 90 == East, 180 == South, 270 == West).
}
\examples{
northangle(0,0,1,1)    # 45
northangle(0,0,-1,-1)  # 225
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}

