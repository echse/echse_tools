\name{skillAsBinary}
\alias{skillAsBinary}
\title{Binary forecast verification measures applied to continous variables}
\usage{
skillAsBinary(x, colObs = "obs", colsSim = c("sim"),
  thresholds = quantile(x[, colObs], probs = c(0.9, 0.95, 0.99)), nMin = 5)
}
\arguments{
  \item{x}{A data frame holding observations and the
  corresponding predictions.}

  \item{colObs}{Name of the column in \code{x} holding the
  observations.}

  \item{colsSim}{A vector of names identifying the columns
  in \code{x} that contain the predictions.}

  \item{thresholds}{A numeric vector defining the
  thresholds for which scores are computed. The default
  uses selected quantiles of the observations.}

  \item{nMin}{The minimum number of observed (POD) or
  predicted (FAR) events to accept the results. If, for a
  particular threshold value, the number of events is less
  than \code{nMin}, the result for this threshold is set to
  \code{NA}. ETS is set to \code{NA} if either of the
  number of observed and predicted events is less than
  \code{nMin}.}
}
\value{
A list with the components 'POD', 'FAR', and 'ETS'. Each of
theses is a data frame with a column 'threshold' followed
by as many columns as present in \code{colsSim}. Column
names are copied from that argument.
}
\description{
The function computes the probability of detection (POD),
the false alarm rate (FAR), and the equitable threat score
(ETS) for a series of observations and one or more
predictions.
}
\note{
For POD and FAR, the range is 0 (worst) to 1 (perfect). The
range of ETS is -1/3 (worst) to 1 (perfect). An ETS of 0
indicates no skill.
}
\examples{
d= data.frame(obs=runif(1000), sim1=runif(1000), sim2=runif(1000))
skillAsBinary(x=d, colsSim=c("sim1","sim2"),  thresholds=c(0.5,0.75,0.9))
}
\author{
David Kneis \email{david.kneis@uni-potsdam.de}
}
\references{
\url{http://www.cawcr.gov.au/projects/verification}
}

