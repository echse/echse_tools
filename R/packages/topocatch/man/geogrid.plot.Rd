\name{geogrid.plot}
\alias{geogrid.plot}
\title{Plot a geogrid object}
\usage{
  geogrid.plot(grid, title = "", randomColors = FALSE,
    len = if (randomColors) 0 else 5, fun = function() { },
    file = NULL, type = "png", replace = FALSE, ...)
}
\arguments{
  \item{grid}{A geogrid object.}

  \item{title}{A title string to appear above the plot's
  legend.}

  \item{randomColors}{If \code{TRUE}, a random color is
  assigned to each unique value. Otherwise, a continuous
  color scale is used.}

  \item{len}{The number of entries in the legend of
  continuous plots. The value is only used if
  \code{randomColors} is \code{TRUE}.}

  \item{fun}{A function containing secondary plot commands
  (like \code{points} or \code{lines} or \code{text}) to
  draw additional objects or annotations on the grid's
  surface. Defaults to a function with an empty body.}

  \item{file}{A file name for the output. If \code{NULL}
  (default), the output appears on the screen only.}

  \item{type}{A file extension string.  Currently 'png' or
  'jpg'.  If \code{file} is \code{NULL}, the value of this
  argument is ignored.}

  \item{replace}{Is is OK to replace an existing output
  file? If \code{file} is \code{NULL}, the value of this
  argument is ignored.}

  \item{...}{Optional arguments to be passed to
  \code{fun}.}
}
\value{
  \code{NULL}
}
\description{
  The geogrid's matrix is converted into a raster image
  that can be either send to a file or the screen.
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}

