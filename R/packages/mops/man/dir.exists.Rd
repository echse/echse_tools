\name{dir.exists}
\alias{dir.exists}
\title{Check existence of a directory}
\usage{
dir.exists(dirname)
}
\arguments{
  \item{dirname}{A character string specifying the single
  path to be checked.}
}
\value{
A logical value.
}
\description{
The function checks whether a \emph{single} directory
exists by evaluating the result of a call to
\code{\link{file.info}}.
}
\author{
David Kneis \email{david.kneis@uni-potsdam.de}
}
\seealso{
See \code{\link{file.exists}} for a similar test on files.
}

