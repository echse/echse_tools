\name{file.move}
\alias{file.move}
\title{Move a single file to another directory}
\usage{
file.move(file, dir, overwrite = FALSE)
}
\arguments{
  \item{file}{The path of the file to be moved (character
  string).}

  \item{dir}{The path of the target directory (character
  string).}

  \item{overwrite}{A logical value. If \code{TRUE},
  existing files in the target directory will silently be
  overwritten. Defaults to \code{FALSE}.}
}
\value{
\code{TRUE} on success and \code{FALSE} otherwise.
}
\description{
The function moves a \emph{single} file to another
directory using the methods \code{\link{file.copy}} and
\code{\link{file.remove}}.
}
\author{
David Kneis \email{david.kneis@uni-potsdam.de}
}
\seealso{
See \code{\link{file.copy}} for copying.
}

