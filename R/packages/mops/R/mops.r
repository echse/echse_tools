#' Utilities to optimize parameters of dynamic simulation models
#'
#' Type \code{help(package="mops")} to inspect the package description.
#'
#' @name mops-package
#' @aliases mops
#' @docType package
{}

.onLoad = function(libname, pkgname) {
  tz= Sys.timezone()
  if (tz != "UTC")
    stop(paste("Current system time zone is '",tz,"' but the '",pkgname,
      "' package requires it to be 'UTC'."))
}

################################################################################

#' Read a time series file
#'
#' The function reads a time series file and returns a data frame,
#' possibly filtered according to supplied rules. The time series is not
#' necessarily regular.
#'
#' @param file Name of the file containing the time series.
#' @param colTime Name of column in \code{file} containing time info.
#' @param colValue Name of column in \code{file} containing the data.
#' @param colsep Column separator used in \code{file}. Defaults to TAB.
#' @param timeConv A function to convert the strings in column
#'   \code{colTime} of \code{file} into values of class \code{POSIXct}.
#'   The default function expects the strings to be in ISO8601 format with date
#'   and time separated by a single blank. Example: '1969-12-31 23:55:00' means
#'   five minutes before the begin of the UNIX epoche. When dealing with longer
#'   time series, it is important to use a time zone without DST, such as 'UTC',
#'   in the \code{tz} argument of the conversion functions. Otherwise, the 
#'   settings for the current locale are likely to be used, which may fail.
#' @param periods.select A data frame with two columns 'begin' and 'end' both
#'   of class \code{POSIXct}. Each record defines a 'period of interest'. The
#'   time series is filtered for data falling into (at least)
#'   one of these periods. Other data are ignored. The limits defined by
#'   'begin' and 'end' are part of the period. If \code{periods} is an empty data
#'   frame (default), no filtering takes place. Note the \code{periods.ignore}
#'   argument.
#' @param periods.ignore Similar to \code{periods.select}, but here, each
#'   record defines a 'period to ignore'. Thus, the time series is filtered for
#'   data \emph{not} falling into \emph{any} of the specified periods. If a
#'   record in \code{periods.ignore} conflicts with a record in
#'   \code{periods.select}, i.e. a 'period to ignore' overlaps with a
#'   'period of interest', data within the overlap are \emph{ignored}.
#' @param nodata A value which indicates missing or corrupt data in the time
#'   series values. The corresponding records are skipped when creating
#'   the result table. Often used values are -9999 or \code{NA} (default).
#'
#' @return A data frame with two columns 'time' and 'value'.
#'
#' @note The function generates an error if the number of records of the
#'   (filtered) time series is zero.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export

read_timeSeries= function(
  file,
  colTime,
  colValue,
  colsep="\t",
  timeConv= function(x) { as.POSIXct(strptime(x,
    "%Y-%m-%d %H:%M:%S",tz="UTC"),tz="UTC") },
  nodata=NA,
  periods.select= data.frame(begin=c(), end=c()),
  periods.ignore= data.frame(begin=c(), end=c())
) {
  # Check args
  if ((length(file) != 1) || (!is.character(file)))
    stop("File name must be a character string.")
  if (!file.exists(file))
    stop(paste("File with time series data '",file,"' not found.",sep=""))
  if ((length(colTime) != 1) || (!is.character(colTime)))
    stop("Name of time column must be a character string.")
  if ((length(colValue) != 1) || (!is.character(colValue)))
    stop("Name of values column must be a character string.")
  if ((length(timeConv) != 1) || (!is.function(timeConv)))
    stop("Need a single function to convert times.")
  if (nrow(periods.select) > 0) {
    if (!all(c("begin","end") %in% names(periods.select)))
      stop("Expecting columns 'begin' and 'end' in data frame 'periods.select'.")
    if ((class(periods.select$begin)[1] != "POSIXct") || (class(periods.select$end)[1] != "POSIXct"))
      stop("Columns 'begin' and 'end' in data frame 'periods.select' must be of class POSIXct.")
    if (any(periods.select$begin > periods.select$end))
      stop("Bad time period specification in data frame 'periods.select'. Check limits.")
  }
  if (nrow(periods.ignore) > 0) {
    if (!all(c("begin","end") %in% names(periods.ignore)))
      stop("Expecting columns 'begin' and 'end' in data frame 'periods.ignore'.")
    if ((class(periods.ignore$begin)[1] != "POSIXct") || (class(periods.ignore$end)[1] != "POSIXct"))
      stop("Columns 'begin' and 'end' in data frame 'periods.ignore' must be of class POSIXct.")
    if (any(periods.ignore$begin > periods.ignore$end))
      stop("Bad time period specification in data frame 'periods.ignore'. Check limits.")
  }
  # Read data
  dat= read.table(file=file, header=TRUE, sep=colsep)
  if (!all(c(colTime, colValue) %in% names(dat))) {
    stop(paste("Missing column(s) in file '",file,"'. Expected times in column '",
      colTime,"' and values in column '",colValue,"'.",sep=""))
  }
  dat= dat[, c(colTime,colValue)]
  dat[,colTime]= timeConv(dat[,colTime])
  if (any(is.na(dat[,colTime])))
    stop(paste("Failed to convert time values read from '",file,"'.",sep=""))
  # Time filter 1 (select filter)
  if (nrow(periods.select) > 0) {
    selected= rep(FALSE, nrow(dat))
    for (i in 1:nrow(periods.select)) {
      selected[(dat[,colTime] >= periods.select$begin[i]) &
        (dat[,colTime] <= periods.select$end[i])]= TRUE
    }
    dat= dat[selected,]
  }
  # Time filter 2 (ignore filter)
  if (nrow(periods.ignore) > 0) {
    selected= rep(TRUE, nrow(dat))
    for (i in 1:nrow(periods.ignore)) {
      selected[(dat[,colTime] >= periods.ignore$begin[i]) &
        (dat[,colTime] <= periods.ignore$end[i])]= FALSE
    }
    dat= dat[selected,]
  }
  # Nodata filter
  if (is.na(nodata)) {
    dat= subset(dat, !is.na(dat[,colValue]))
  } else {
    dat= subset(dat, dat[,colValue] != nodata)
  }
  # Check length
  if (nrow(dat) == 0)
    stop(paste("Filtered time series read from '",file,"' has zero records.",sep=""))
  # Set colnames
  names(dat)[which(names(dat) == colTime)]= "time"
  names(dat)[which(names(dat) == colValue)]= "value"
  # Return data frame
  return(dat)
}

################################################################################

#' Run an external model
#'
#' The function runs a model (being an executable file) with the supplied set
#' of command line arguments and captures the return code.
#'
#' @param model_path Path of the model executable. Note that the executable is
#'   executed in R's current working directory, i.e. the one returned by
#'   \code{\link{getwd}}. This is true even if the executable resides
#'   somewhere else.
#' @param model_args A \emph{named} vector of command line arguments to be
#'   passed to the model. Note that, for named elements, name and value are
#'   concatenated by '=', automatically. Moreover, all values are enclosed in
#'   shell quotes (using \code{\link{shQuote}}). For example, a vector
#'   \code{c(a='1', '-x', b='2')} is expanded to the string 'a='1' '-x' b='2''.
#' @param stopIfNonZero If \code{TRUE} the function calls \code{stop} if the
#'   model executable returns a non-zero exit code. If \code{FALSE}, only a
#'   warning is generated.
#'
#' @return The return code of the model executable.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export

run_model= function(model_path, model_args, stopIfNonZero=FALSE) {
  if (is.null(names(model_args))) {
    sep= rep("",length(model_args))
  } else {
    sep= rep("=",length(model_args))
    sep[names(model_args) == ""]= ""
  }
  cmd= paste(model_path," ",paste(names(model_args),sep,shQuote(model_args),
    sep="",collapse=" "),sep="")
  status= system(cmd, intern=FALSE, ignore.stderr=FALSE, wait=TRUE)
  if (status != 0) {
    if (stopIfNonZero) {
      stop(paste("Model run failed with code ",status,". The command line was: '",cmd,"'.",sep=""))
    } else {
      warning(paste("Model run failed with code ",status,". The command line was: '",cmd,"'.",sep=""))
    }
  }
  return(status)
}

################################################################################

#' Compute the error of a dynamic model simulation for a given set of parameters
#'
#' The function calls a dynamic simulation model (an external executable) with a
#' given set of model parameters. It then computes measure(s) of the model error
#' by comparing simulated and observed time series.
#' The function takes the vector of model parameters as the first argument,
#' hence it can be used as an objective function in R's optimization method
#' \code{\link{optim}}.
#'
#' @param parameters A \emph{named} vector holding the values of the parameters
#'    to be used in the model call. Note that necessary changes to the model's input
#'    file(s) are \emph{NOT} performed automagically. Any of those changes need to be
#'    done by the user-specified function \code{func_first} (see below).
#' @param func_first An external function to be called at the very beginning 
#'   of the current function.  It is called as \code{f(parameters, moreArgs_first)}
#'   where \code{parameters} is the parameters vector and \code{moreArgs_first}
#'   holds any additional arguments. 
#'   Typically, \code{func_first} is responsible for changing the model's input
#'   files to reflect the current parameter set contained in \code{parameters}. This
#'   can be achieved, for example, by calls to \code{\link{update_template}} or
#'   similar methods (see examples). Furthermore, \code{func_first} often performs an
#'   initial clean-up if \code{modelError_multiDim} is called repeatedly under
#'   control of another method (such as \code{\link{optim}}, for example).
#'   The function should return an integer code of 0 to signal successful execution
#'   or a non-zero code to signal an error. 
#' @param moreArgs_first Last argument to be passed to \code{func_first}.
#'   Use a list if the function requires multiple values as inputs.
#' @param func_final An external function to be called at the very end of the
#'   current function. It is called as \code{f(parameters, gof, moreArgs_final)}
#'   where \code{gof} is the return value of \code{gof_function} (see below).
#'   Further arbitrary arguments can be passed via \code{moreArgs_final}.
#'   If \code{modelError_multiDim} is called repeatedly under the control of
#'   another method, \code{func_final} is typically responsible for saving the
#'   model's output and/or recording the godness-of-fit for the tested parameter
#'   sets. See argument \code{func_first} for possible return values.
#' @param moreArgs_final Last argument to be passed to \code{func_final}.
#'   Use a list if multiple values need to be passed.
#' @param observed A data frame holding one or more time series of observations.
#'   The time column specified as \code{obs_colTime} must be of type
#'   \code{POSIXct}, other columns must be \code{numeric}.
#'   Any Missing values must be marked with \code{NA}.
#' @param obs_colTime Column in table \code{observed} holding time information.
#' @param sim_files A named vector of file names specifying the simulated time
#'   series to be compared with observation data. Element names must correspond
#'   to column names in \code{observed}.
#' @param sim_colsTime Vector of the same length as \code{sim_files} specifying
#'   the names of the columns containing time info.
#' @param sim_colsValue Vector of the same length as \code{sim_files} specifying
#'   the names of the columns containing the data.
#' @param sim_colsep Column separator used in \code{sim_file}. Defaults to TAB.
#' @param sim_timeConv A function to convert the strings in column
#'   \code{sim_colTime} of \code{sim_file}. See argument \code{timeConv} of
#'   the \code{\link{read_timeSeries}} method for details.
#' @param periods.select A data frame with a set of temporal filter rules being applied
#'   to the time series. See the \code{periods.select} argument of the
#'   \code{\link{read_timeSeries}} method for details.
#' @param periods.ignore A data frame with a set of temporal filter rules being applied
#'   to the time series. See the \code{periods.ignore} argument of the
#'   \code{\link{read_timeSeries}} method for details.
#' @param gof_function An objective function of the form 'f(obs,sim)' which returns the
#'   goodness-of-fit based on the vectors of observed (obs) and simulated data
#'   (sim), respectively. In the context of optimization, it should return a
#'   scalar result (see the documentation of the \code{\link{optim}} method). If
#'   the result is not a scalar, it should be a vector, preferrably with names.
#' @param unlist If \code{TRUE}, the \code{\link{unlist}} method is applied to
#'   the return value. Default is \code{FALSE}.
#' @param stopIfNoData If \code{TRUE}, an error is generate if no corresponding
#'   observed and simulated data are found. If \code{FALSE}, the function
#'   \code{gof_function} is called as usual but both arguments are vectors of
#'   zero length, thus, \code{gof_function} must handle this special case.
#' @inheritParams run_model
#'
#' @return By default, a list of the same length (and with the same names) as vector \code{files_sim}.
#'   Each list element holds the return value of \code{gof_function} applied to
#'   the time series of simulated and observed values. If the \code{unlist}
#'   argument is \code{TRUE}, the list is converted into a vector (using the
#'   \code{\link{unlist}} method).
#'
#' @note
#'   For 1-dimensional optimization using R's \code{\link{optimize}}
#'   method, use \code{\link{modelError_oneDim}} which is a convenient wrapper
#'   for \code{modelError_multiDim}.
#'
#' @seealso \code{\link{modelError_oneDim}} for the 1-dimensional case
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Example for 'func_first' that updates parameter values in a single file
#' updatePars_byTemplate= function(p, args) {
#'   stopifnot(all(c("fileTemplate","fileResult","charOpen","charClose")
#'     %in% names(args)))
#'   n= update_template(file_template=args$fileTemplate, vect_updates=p,
#'     char_open=args$charOpen, char_close=args$charClose,
#'     file_result=args$fileResult, overwrite=TRUE)
#'   if (n == 0)
#'     stop(paste("Updating of '",args$fileTemplate,"' failed.",sep=""))
#' }
#' 
#' # Example for 'func_first' that updates parameter values in multiple files
#' updatePars_byTemplates= function(p, args) {
#'   stopifnot(all(c("filesTable","charOpen","charClose") %in% names(args)))
#'   stopifnot(is.data.frame(args$filesTable),
#'    all(c("fileTemplate","fileResult") %in% names(args$filesTable)))
#'   for (i in 1:nrow(args$filesTable) {
#'     n= update_template(file_template=args$filesTable$fileTemplate[i],
#'       vect_updates=p, char_open=args$charOpen, char_close=args$charClose,
#'       file_result=args$filesTable$fileResult[i], overwrite=TRUE)
#'     if (n == 0)
#'       stop(paste("Updating of '",args$filesTable$fileTemplate[i],
#'         "' failed.",sep=""))
#'   }
#' }
#'
#' }

modelError_multiDim = function(
  parameters,
  model_path,
  model_args,
  func_first,
  moreArgs_first=NULL,
  func_final,
  moreArgs_final=NULL,
  observed,
  obs_colTime,
  sim_files,
  sim_colsTime,
  sim_colsValue,
  sim_colsep="\t",
  sim_timeConv= function(x) { as.POSIXct(strptime(x,
    "%Y-%m-%d %H:%M:%S",tz="UTC"),tz="UTC") },
  periods.select= data.frame(begin=c(), end=c()),
  periods.ignore= data.frame(begin=c(), end=c()),
  gof_function= function(obs,sim) { mean((obs-sim)^2) },
  unlist=FALSE,
  stopIfNoData= TRUE
) {
  #.............................................................................
  # Check args
  tryCatch({
    stopifnot(length(parameters)>0, (!is.null(names(parameters))), all(names(parameters) != ""))
    stopifnot(is.function(func_first))
    stopifnot(is.function(func_final))
    stopifnot(is.data.frame(observed), ncol(observed)>=2, nrow(observed)>0,
      class(observed[,obs_colTime])[1]=="POSIXct")
    stopifnot(length(sim_files) >= 1,
     !any(is.null(names(sim_files))), all(names(sim_files) != ""),
     length(sim_colsTime) == length(sim_files),
     length(sim_colsValue) == length(sim_files))
  }, error= function(e) {
    stop(paste("Error in arguments. Details: ",e,sep=""))
  })
  #.............................................................................
  # Call init function
  tryCatch({
    stat= func_first(parameters, moreArgs_first)
    if ((!is.numeric(stat)) || (stat != 0))
      stop("'func_first' returned a non-zero exit code.")
  }, error= function(e) {
    stop(paste("Error in initialization function 'func_first'. Details: ",e,sep=""))
  })
  #.............................................................................
  # Run model
  run_model(model_path, model_args, stopIfNonZero=TRUE)
  #.............................................................................
  # Import model outputs
  tryCatch({

    # Read all simulated time series
    for (i in 1:length(sim_files)) {
      tmp= read_timeSeries(file=sim_files[i], colTime=sim_colsTime[i], colValue=sim_colsValue[i],
        colsep=sim_colsep, timeConv= sim_timeConv, nodata=NA,
        periods.select=periods.select, periods.ignore=periods.ignore)
      if (i == 1) {
        sim= tmp[,c("time", "value")]
      } else {
        if (!identical(tmp$time, sim$time))
          stop(paste("Time column in '",sim_files[i],"' differs from time column in '",sim_files[1],"'.",sep=""))
        sim= cbind(sim, tmp[,"value"])
      }
      names(sim)[ncol(sim)]= names(sim_files)[i]
    }
    rm(tmp)
  }, error= function(e) {
    stop(paste("Failed to collect simulated data. Details: ",e,sep=""))
  })

  #.............................................................................
  # Apply the error function to the data

  # Make sure that there is a matching column with observations for every simulated time series
  if (!all(names(sim)[2:ncol(sim)] %in% names(observed)[which(names(observed) != obs_colTime)]))
    stop("Element names of vector 'sim_files' must match column names in the table of observations.")

  # Truncate the time axes of obs and sim to the range with possible matches
  tmin= max(c(min(observed[,obs_colTime]), min(sim$time)))
  tmax= min(c(max(observed[,obs_colTime]), max(sim$time)))
  observed= subset(observed, ((observed[,obs_colTime] >= tmin) & (observed[,obs_colTime] <= tmax)))
  sim= subset(sim, ((sim$time >= tmin) & (sim$time <= tmax)))

  tryCatch({
    gof= list()
    for (i in 2:ncol(sim)) {
      tmp= merge(x=observed[,c(obs_colTime,names(sim)[i])], y=sim[,c(1,i)],
        by.x=obs_colTime, by.y="time", all=FALSE, suffixes= c(".obs",".sim"))
      colObs= paste(names(sim)[i],"obs",sep=".")
      colSim= paste(names(sim)[i],"sim",sep=".")
      tmp= subset(tmp, (!is.na(tmp[,colObs]) & !is.na(tmp[,colSim])))
      if ((nrow(tmp) == 0) && (stopIfNoData))
        stop(paste("No corresponding observations for simulated data with name '",names(sim)[i],"'.",sep=""))
      gof[[i-1]]= gof_function(obs=tmp[,colObs], sim=tmp[,colSim])
      names(gof)[i-1]= names(sim)[i]
    }
    rm(tmp, colObs, colSim)
  }, error= function(e) {
    stop(paste("Failed to compute model error. Details: ",e,sep=""))
  })

  #.............................................................................
  # Call finalize function
  tryCatch({
    stat= func_final(parameters, gof, moreArgs_final)
    if ((!is.numeric(stat)) || (stat != 0))
      stop("'func_final' returned a non-zero exit code.")
  }, error= function(e) {
    stop(paste("Error in finalize function 'func_final'. Details: ",e,sep=""))
  })
  #.............................................................................
  # Return the computed error
  if (unlist) {
    return(unlist(gof))
  } else {
    return(gof)
  }
}


################################################################################

#' Compute the error of a dynamic model simulation for a given parameter value
#'
#' This is a wrapper for \code{\link{modelError_multiDim}} for the purpose of
#' 1-dimensional optimization using R's \code{\link{optimize}} method. Instead
#' of a named vector of parameters, it takes the value and name of a single
#' parameter. See \code{\link{modelError_multiDim}} for more details.
#'
#' @param p.value Value of the single model parameter which is to be updated
#'   before calling the model.
#' @param p.name Name of the parameter in \code{p.value}.
#' @inheritParams modelError_multiDim
#'
#' @return See return value of \code{\link{modelError_multiDim}}.
#'
#' @seealso Use \code{\link{modelError_multiDim}} for cases where more than one
#'   model parameter is to be varied as in multi-dimensional optimization,
#'   for example.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export

modelError_oneDim= function(
    p.value,
    p.name,
    model_path,
    model_args,
    func_first,
    moreArgs_first=NULL,
    func_final,
    moreArgs_final=NULL,
    observed,
    sim_files,
    sim_colsTime,
    sim_colsValue,
    sim_colsep="\t",
    sim_timeConv= function(x) { as.POSIXct(strptime(x,
      "%Y-%m-%d %H:%M:%S",tz="UTC"),tz="UTC") },
    periods.select= data.frame(begin=c(), end=c()),
    periods.ignore= data.frame(begin=c(), end=c()),
    gof_function= function(obs,sim) { mean((obs-sim)^2) },
    unlist=FALSE,
    stopIfNoData=TRUE
  ) {
    p= c(p.value)
    names(p)= p.name
  modelError_multiDim(parameters=p,
    model_path=model_path, model_args=model_args,
    func_first=func_first, moreArgs_first=moreArgs_first,
    func_final=func_final, moreArgs_final=moreArgs_final,
    observed=observed,
    sim_files=sim_files, sim_colsTime=sim_colsTime, sim_colsValue=sim_colsValue,
    sim_colsep=sim_colsep, sim_timeConv=sim_timeConv,
    periods.select=periods.select, periods.ignore=periods.ignore,
    gof_function=gof_function,
    unlist= unlist,
    stopIfNoData= stopIfNoData
  )
}

################################################################################

#' Evaluate an objective function for a number of pre-defined parameter sets 
#'
#' The method evaluates an objective function for a fixed number of parameter
#' sets. The parameter sets are generating by building all possible combinations
#' of the test-values supplied for the individual parameters.
#'
#' @param p A list defining the parameter sets to be tested. Each
#'   element of the list should be a numeric vector, representing the test values
#'   of a particular parameter. List elements must be named. The names are passed
#'   as the parameters' names to \code{f}.
#' @param f Objective function with interface \code{f(par,...)}.
#'   The 1st argument of this function must be a numeric vector holding the
#'   parameters. When \code{f} is called, the names of the elements in \code{par}
#'   are copied from the names of the elements in argument \code{p}. The return
#'   value of \code{f} should be a \emph{named} vector (which can be of length 1).
#' @param ... Additional arguments to be passed to the objective function \code{f}. 
#'
#' @return A list with two components:
#'   \item{params}{A data frame holding all tested parameter sets. The column
#'     names are identical to the names of \code{p}.}
#'   \item{errors}{A list whose length equals the number of tested parameter sets.
#'     The list element with index k holds the return value of \code{f} for
#'     the parameter set stored in the k-th row of the returned data
#'     frame \code{params}. The \code{\link{lapply}} method is suitable for
#'     extracting values from the list.}
#'
#' @seealso Use \code{\link{mcs_run}} which uses random parameter sets instead of
#'   pre-defined ones.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#'
#' @examples
#' # Approximate the parameters of a parabola y=ax^2+b by simple grid search
#' x.obs=c(0,1,2)
#' y.obs=c(0,1,4)
#' objfun= function(p, x.obs, y.obs) {
#'   c(sse= sum(((p["a"]*x.obs^2+p["b"]) - y.obs)^2))
#' }
#' p= list(a=seq(-1,1,0.1), b=seq(-1,1,0.1))
#' res= gridSearch(p=p, f=objfun, x.obs=x.obs, y.obs=y.obs)
#' sse= unlist(lapply(X=res$errors, FUN=function(z) { z[["sse"]] }))
#' print(res$params[which.min(sse),])

gridSearch= function(p, f, ...) {
  # check args
  if (!is.function(f))
    stop("Argument 'f' must be a function.")
  if ((!is.list(p)) || is.null(names(p)) || any(names(p) == ""))
    stop("Argument 'p' must be a list where all elements are named.")
  if (any(duplicated(names(p))))
    stop("Found non-unique element names in list supplied as argument 'p'.")
  for (i in 1:length(p)) {
    if (!is.numeric(p[[i]]))
      stop(paste("Element '",names(p)[i],"' of argument 'p' is not numeric.",sep=""))
    if (length(p[[i]]) < 1)
      stop(paste("Element '",names(p)[i],"' of argument 'p' is empty.",sep=""))
  }
  # generate data frame with all parameters sets (each row is 1 set)
  sets= expand.grid(p, KEEP.OUT.ATTRS=FALSE)
  # evaluate objective function for all sets
  errs= list()
  for (i in 1:nrow(sets)) {
    pars= as.numeric(unlist(sets[i,]))
    names(pars)= names(sets)
    errs[[i]]= f(pars, ...)
  }
  return(list(params=sets, errors=errs))
}

