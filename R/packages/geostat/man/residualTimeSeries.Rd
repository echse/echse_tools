\name{residualTimeSeries}
\alias{residualTimeSeries}
\title{Computes time series of residuals assuming a linear model}
\usage{
  residualTimeSeries(file_obs, file_loc, file_resid,
    file_coeff, colsep = "\t", obs_colTime = "datetime",
    loc_colId = "id", loc_colPred = "z", average = FALSE,
    r2min = 0.36, nsignif = 3, overwrite = FALSE)
}
\arguments{
  \item{file_obs}{Input file containing the time series of
  observations. Column names are expected in the first row.
  The column names are interpreted as location IDs except
  for the one specified by \code{obs_colTime} (header of
  the time column). The strings in the time column can be
  in any convenient format, for example '2000-01-01
  12:30:00'. The time column is copied to the output files
  'as is'.}

  \item{file_loc}{Input file containing a locations'
  attribute. Must be a table with two named columns at
  least, one holding the locations' IDs and the other one
  containing the value of a numerical attribute (typically
  elevation).  For each location in \code{file_obs}, there
  must be one matching record.}

  \item{file_resid}{Output file holding the time series of
  residuals. The file's layout is identical with that of
  the input file \code{file_obs}.}

  \item{file_coeff}{Output file holding the time series of
  coefficients of the linear model. The first column in
  this file is identical with the time column from
  \code{file_obs}. The columns 'intercept' and 'slope' hold
  the model coefficients and the column 'r2' list the
  values of R-squared.}

  \item{colsep}{Field separator used in all input and
  output files.  Defaults to TAB.}

  \item{obs_colTime}{Name of the column in \code{file_obs}
  containing time info instead of numerical data. Defaults
  to 'datetime'.}

  \item{loc_colId}{Name of the column in \code{file_loc}
  holding the locations' IDs. Defaults to 'id'.}

  \item{loc_colPred}{Name of the column in \code{file_loc}
  holding the values of the predictor variable. Defaults to
  'z'.}

  \item{average}{Logical argument (default \code{FALSE}).
  If \code{FALSE}, regression coefficients are computed for
  individual records, i.e. for every single time step. If
  \code{TRUE}, the data are averaged over time before the
  coefficients of the linear model are estimated. These
  parameters are then applied to \emph{all} time steps.}

  \item{r2min}{Minimum value of R-squared for the linear
  model relating the observed variable to the predictor
  variable. The linear model is only used if the value of
  R-squared for a particular time step is >= \code{r2min}.
  In case of a lower value of R-squared, slope and
  intercept of the linear model are set to zero and the
  residuals become identical with the original data.  The
  default for \code{r2min} is 0.36 which is equivalent to a
  coefficient of correlation of 0.6.}

  \item{nsignif}{The number of \emph{significant} digits to
  be used when printing the residuals to \code{file_resid}.
  Defaults to 3.}

  \item{overwrite}{Defaults to FALSE. If TRUE, existing
  output files are replaced without a warning.}
}
\value{
  \code{NULL}
}
\description{
  The function takes a multi-location time series of an
  observed variable \eqn{y} as input. Typically, \eqn{y} is
  a meteorological variable. In addition, the function
  requires information on a static predictor variable
  \eqn{z}, whose value is known at the locations where
  \eqn{y} has been measured. Assuming the existance of a
  \emph{linear} correlation \eqn{y= a * z + b}, the
  function converts the original data (\eqn{y}) into
  (additive) residuals \eqn{r = y - (a * z + b)}. The
  coefficients \eqn{a} and \eqn{b} are estimated either for
  every single time step or time-averaged data. The
  function returns both a time series of the residuals
  \eqn{r} as well as a time series of the linear model's
  coefficients \eqn{a} and \eqn{b}.
}
\note{
  When using the output of this function for spatial
  residual interpolation, one needs to take care of
  possible undesired extrapolation effects. Problems are
  likely to occur if the interpolated variable is (1)
  physically limited to a certain range and (2) the value
  of the predictor variable at a target location is far
  outside the predictor variable's range at the source
  locations. Example: Rainfall is interpolated for a
  low-elevation site using rain gage data from higher
  sites. If the correlation between elevation and rainfall
  is strong but negative, the estimated rainfall for the
  low-elevation site may be negative.
}
\examples{
# Test input data
print("__________ Obs. data __________")
obs= data.frame(
  datetime=ISOdatetime(2000,1,1,0,0,0) + seq(from=0, to=2, by=1) * 3600,
  stationA=c(1,2,8), stationB=c(2,3,6), stationC=c(6,2,1), stationD=c(4,3,3))
print(obs)
write.table(obs, file=fObs<-tempfile(), col.names=TRUE,row.names=FALSE,sep="\\t")
print("__________ Locations __________")
loc= data.frame(
  id= c("stationA","stationB","stationC","stationD"),
  elevation= c(2, 400, 1800, 900))
print(loc)
write.table(loc, file=fLoc<-tempfile(), col.names=TRUE,row.names=FALSE,sep="\\t")
# Convert
residualTimeSeries(
  file_obs= fObs, file_loc= fLoc,
  file_resid=fResid<-tempfile(), file_coeff=fCoeff<-tempfile(),
  colsep="\\t",
  obs_colTime="datetime", loc_colId="id", loc_colPred="elevation",
  nsignif=4, overwrite=FALSE
)
# Result of conversion
print("__________ Residuals __________")
resid= read.table(file=fResid, header=TRUE, sep="\\t")
print(resid)
print("__________ Model fit __________")
coeff= read.table(file=fCoeff, header=TRUE, sep="\\t")
print(coeff)
# Back conversion of residuals
print("__________ Back conv. __________")
back= resid
for (i in 2:ncol(resid))
  back[,i]= resid[,i] +
    (coeff$slope * loc$elevation[match(names(obs)[i],loc$id)] + coeff$intercept)
print(back)
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}

