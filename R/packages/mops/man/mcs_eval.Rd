\name{mcs_eval}
\alias{mcs_eval}
\title{Evaluation of Monte-Carlo simulation output}
\usage{
mcs_eval(outdir_mcs, obs_files, obs_colsTime, obs_colsValue,
  obs_colsep = "\\t", obs_timeConv = function(x) {    
  as.POSIXct(strptime(x, "\%Y-\%m-\%d \%H:\%M:\%S", tz = "UTC"), tz = "UTC") },
  sim_file, sim_colTime, sim_colValue, sim_colsep = "\\t",
  sim_timeConv = function(x) {     as.POSIXct(strptime(x,
  "\%Y-\%m-\%d \%H:\%M:\%S", tz = "UTC"), tz = "UTC") }, obs_nodata,
  gof_function = function(obs, sim) {     c(bias = mean(sim - obs), rmse =
  sqrt(mean((sim - obs)^2))) })
}
\arguments{
  \item{outdir_mcs}{Name of the directory created by
  \code{\link{mcs_run}} containing the output of a
  Monte-Carlo simulation.}

  \item{obs_files}{A \emph{named} vector of file names.
  Each of the files contains a time series of observations
  against which simulation results are compared.}

  \item{obs_colsTime}{A vector of the same length as
  \code{obs_files}, specifying the name(s) of the colum(s)
  holding time information.}

  \item{obs_colsValue}{A vector of the same length as
  \code{obs_files}, specifying the name(s) of the colum(s)
  holding the observation data.}

  \item{sim_file}{\emph{Base name} of the file containing
  the simulated time series. The full file path is
  assembled automatically.}

  \item{sim_colsep}{Column separator used in
  \code{sim_file}. Defaults to TAB.}

  \item{sim_timeConv}{A function to convert the strings in
  column \code{sim_colTime} of \code{sim_file}. See
  argument \code{timeConv} of the
  \code{\link{read_timeSeries}} method for details.}

  \item{gof_function}{An objective function of the form
  'f(obs,sim)' which returns the goodness-of-fit based on
  the vectors of observed (obs) and simulated data (sim),
  respectively. In the context of optimization, it should
  return a scalar result (see the documentation of the
  \code{\link{optim}} method). If the result is not a
  scalar, it should be a vector, preferrably with names.}

  \item{obs_colsep}{Column separator used in
  \code{obs_file}. Defaults to TAB.}

  \item{obs_timeConv}{A function to convert the strings in
  column \code{obs_colTime} of \code{obs_file}. See
  argument \code{timeConv} of the
  \code{\link{read_timeSeries}} method for details.}

  \item{sim_colTime}{Name of column in \code{sim_file}
  containing time info.}

  \item{sim_colValue}{Name of column in \code{sim_file}
  containing the data.}

  \item{obs_nodata}{A value which indicates missing or
  corrupt data in the time series of observations. Often
  used values are -9999 or \code{NA}.}
}
\value{
The function returns an object of class 'mcs_eval'. This is
a list with the following components: \item{tbl}{A data
frame that lists the tested parameter sets together with
the corresponding values of the objective function(s).  The
columns 'err.obs', 'err.len', 'err.fun', and 'err.val'
specify the name of the observation data set, the number of
observations, the name of the objective function, and the
value of the objective function, respectively.}
\item{par}{A vector holding the names of the parameters as
appearing in the column headers of \code{tbl}.}
}
\description{
The function analyzes the output created by
\code{\link{mcs_run}}. It computes goodness-of-fit measures
for all tested parameter sets.
}
\note{
Multiple observation files may be supplied in the
\code{obs_files} argument. This makes it possible to
analyze the relation between the model parameters and the
simulation error for different sub-sets of observations.
This is often useful if the method of Monte-Carlo
simulation is applied in order to identify a set of optimum
model parameters.  For example, in hydrological model
calibration, \code{obs_files} could point to three files:
(1) a time series of all observations, (2) a time series of
observations during winter time and (3) another time series
of observations during low-flow periods. This often
facilitates the identification of those model parameters
which are sensitive in certain periods of time only (e.g.
during snow cover or low flow).
}
\author{
David Kneis \email{david.kneis@uni-potsdam.de}
}

