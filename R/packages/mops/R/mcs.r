
################################################################################
################################################################################
################################################################################

#' Sample generation for Monte-Carlo simulation
#'
#' The function generates a number of random parameter sets using a latin
#' hypercube algorithm from package 'lhs'. The underlying distribution is
#' uniform.
#'
#' @param ranges_table Data frame with 3 columns 'parameter', 'min', 'max'. Each
#'   record defines the name and the sampling range for a particular parameter.
#' @param n The number of parameter sets to be generated.
#'
#' @return A data frame with \code{n} rows and as many columns as there are
#'   records in \code{ranges_table}. The column are named after the parameters.
#'
#' @note If the generated random sets are to be used by a simulation model, the
#'   sampling ranges must be chosen carefully.
#'   In models with many parameters, a particular
#'   value of one parameter may be incompatible with a particular value of
#'   another parameter. For example, if the sampling range
#'   for a parameter 'minimumCapacity' was 2...10, it would not make sense to
#'   use a sampling range 5...20 for a related parameter 'maximumCapacity'. If
#'   the ranges were defined as above, it is likely that 'minimumCapacity' > 
#'   'maximumCapacity' in some sets. Those cases typically lead to
#'   implausible model outputs or a crash of the model due to numeric errors.
#'   In addition, one should be aware of the fact that the initial value(s) of
#'   a model's state variable(s) must be compatible with all possible random
#'   parameter sets. For example, an intial value of 0 for a state variable
#'   'temperature' would be critical if the sampling range for a related
#'   parameter 'minimumTemperature' was 5...50.
#'   
#' @seealso Use \code{\link{mcs_run}} to run a Monte-Carlo simulation using
#'   the generated sets.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#'
#' @examples
#' mcs_sample(data.frame(parameter=c("a","b"),min=c(1,5),max=c(2,10)), n=10)

mcs_sample= function(
  ranges_table,
  n=100
) {
  #.............................................................................
  # Check args
  if (!is.data.frame(ranges_table))
    stop("Table of parameter ranges must be a data frame.")
  if (nrow(ranges_table) == 0)
    stop("Empty table of parameter ranges.")
  if (!all(c("parameter","min","max") %in% names(ranges_table)))
    stop("Missing column(s) in table of parameter ranges.")
  if (length(unique(ranges_table$parameter)) != nrow(ranges_table))
    stop("Non-unique parameter names in table of parameter ranges.")
  ranges_table$parameter= as.character(ranges_table$parameter)
  if (!identical(ranges_table$parameter, make.names(ranges_table$parameter)))  
    stop("Parameter names in table of parameter ranges must be valid R names.")
  if (any(ranges_table$min > ranges_table$max))
    stop("Bad specification of min/max value(s) in table of parameter ranges.")
  if ((!is.numeric(n)) || (length(n) != 1) || (n < 1))
    stop("Invalid number of realizations.")
  #.............................................................................
  # Draw LHS sample
#  if (!library("lhs", character.only=TRUE, quietly=TRUE, logical.return=TRUE))
#    stop("Required package 'lhs' cannot be loaded.")
  parameters= lhs::improvedLHS(n=n, k=nrow(ranges_table))
  #.............................................................................
  # Convert to the actual parameters
  for (i in 1:nrow(parameters)) {
    parameters[i,]= qunif(p=parameters[i,], min=ranges_table$min, max=ranges_table$max)
  }
  parameters= as.data.frame(parameters)
  names(parameters)= ranges_table$parameter
  return(parameters)
}

################################################################################
################################################################################
################################################################################

#' Monte-Carlo simulation for a dynamic simulation model
#'
#' The function runs an external simulation model for each set in a list of
#' parameter sets. The corresponding model output is saved
#' to a directory with a unique time stamp. This time stamp is written to a log
#' file along with the corresponding values of the parameters (see details below).
#' This log file provides the link between parameter values and simulation results
#' which is needed for later analysis.
#'
#' @param sets A data frame of parameter sets. Use \code{\link{mcs_sample}}
#'   to generate a suitable data frame.
#' @param updating_table A data frame with 2 named columns 'file_template' and 'file_result'.
#    Each file in the column 'file_template' is processed with
#    the \code{\link{update_template}} function to create the corresponding
#    result file by replacing placeholders by the appropriate parameter values.
#    For details, refer to the arguments of the \code{\link{update_template}}
#    function with identical names. If all updateable parameter values are read
#    from a single file, the number of rows in the table is 1 (a typical case).
#' @param placeholder_openChar See arguments of \code{\link{update_template}}.
#' @param placeholder_closeChar See arguments of \code{\link{update_template}}.
#' @param outdir_model The directory where the model puts \emph{all} of its
#'   output files. Must be an existing, empty directory.
#' @param outdir_mcs The directory where the results if the Monte-Carlo
#'   simulation are to be collected. Must be an existing, empty directory.
#' @param silent If \code{TRUE}, info on the progress of the Monte-Carlo
#'   simulation is printed. Default is \code{FALSE}.
#'
#' @inheritParams modelError_multiDim
#' @inheritParams run_model
#'
#' @return \code{NULL}
#'
#' @note The log file is created in the directory supplied as \code{outdir_mcs}
#'   and its fixed name is 'mcs.log'. Its contains a table with the time stamp
#'   in the first column. The remaining \code{nrow(ranges_table)} columns hold
#'   the values of the parameters.
#'
#'   In many cases, care must be taken when chosing the parameter ranges (see
#'   argument \code{ranges_table}). In models with many parameters, a particular
#'   value of one parameter may be incompatible with a particular value of
#'   another parameter. This must be taken into account when defining the 
#'   sampling ranges for the two parameters. For example, if the sampling range
#'   for a parameter 'minimumCapacity' was 2...10, it would not make sense to
#'   use a sampling range 5...20 for a related parameter 'maximumCapacity'. If
#'   the ranges were defined as above, it is likely that 'minimumCapacity' > 
#'   'maximumCapacity' in some parameter sets. Those cases typically lead to
#'   implausible model outputs or a crash of the model due to numeric errors.
#'   In addition, one should be aware of the fact that the initial value(s) of
#'   the model's state variable(s) must be compatible with all possible random
#'   parameter sets. For example, an intial value of 0 for a state variable
#'   'temperature' would be critical if the sampling range for a related
#'   parameter 'minimumTemperature' was 5...50.
#'   
#' @seealso Use \code{\link{mcs_eval}} to analyze the output.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export

mcs_run = function(
  sets,
  updating_table,
  placeholder_openChar="{",
  placeholder_closeChar="}",
  model_path,
  model_args,
  outdir_model,
  outdir_mcs,
  silent=FALSE
) {
  # Check parameter table
  if (!is.data.frame(sets))
    stop("Table of parameter sets must be a data frame.")
  if (nrow(sets) == 0)
    stop("Empty table of parameter sets.")
  if (any(duplicated(names(sets))))
    stop("Non-unique names in table of parameter sets.")
  # Check updating table
  if (!is.data.frame(updating_table))
    stop("Table of updating rules must be a data frame.")
  if (nrow(updating_table) == 0)
    stop("Empty table of updating rules.")
  if (!all(c("file_template","file_result") %in% names(updating_table)))
    stop("Missing column(s) in table of updating rules.")
  # Check model settings
  if ((!is.character(model_path)) && (length(model_path) != 1))
    stop("Argument 'model_path' must be a character string.")
  if ((!is.character(model_args)) && (length(model_args) != 1))
    stop("Argument 'model_args' must be a character string.")
  # Check directories
  if ((!is.character(outdir_model)) && (length(outdir_model) != 1))
    stop("Argument 'outdir_model' must be a character string.")
  if (!dir.exists(outdir_model))
    stop("Directory '",outdir_model,"' supplied as 'outdir_model' does not exist.")
  if (length(list.files(path=outdir_model)) > 0)
    stop("Directory '",outdir_model,"' supplied as 'outdir_model' is not empty.")
  if ((!is.character(outdir_mcs)) && (length(outdir_mcs) != 1))
    stop("Argument 'outdir_mcs' must be a character string.")
  if (!dir.exists(outdir_mcs))
    stop("Directory '",outdir_mcs,"' supplied as 'outdir_mcs' does not exist.")  
  if (length(list.files(path=outdir_mcs)) > 0)
    stop("Directory '",outdir_mcs,"' supplied as 'outdir_mcs' is not empty.")
  # More checks
  if ((length(silent) != 1) || (!is.character(silent)))
    silent=FALSE
  #.............................................................................
  # Initialize log file
  logfile= paste(outdir_mcs,"mcs.log",sep="/")
  write(x=paste("timestamp","\t",paste(names(sets),collapse="\t"),sep=""),
    file=logfile, ncolumns=1, append=FALSE)
  #.............................................................................
  # Run MCS
  for (i in 1:nrow(sets)) {
    if (!silent) print(paste("Set ",i," of ",nrow(sets),sep=""))
    parameters= as.numeric(unlist(sets[i,]))
    names(parameters)= names(sets)
    # Update parameter values in template files
    updating_table$file_template= as.character(updating_table$file_template)
    updating_table$file_result= as.character(updating_table$file_result)
    tryCatch({
      for (i in 1:nrow(updating_table)) {
        n= update_template(
          file_template=updating_table$file_template[i],
          vect_updates=parameters,
          char_open=placeholder_openChar, char_close=placeholder_closeChar,
          file_result=updating_table$file_result[i], overwrite=TRUE)
        if (n == 0) stop(paste("Updating of file '",updating_table$file_template[i],
          "' failed. No matching placeholders found.",sep=""))
      }
    }, error= function(e) {
      stop(paste("Updating of parameter file(s) failed. Details: ",e,sep=""))
    })
    # Run model
    run_model(model_path, model_args, stopIfNonZero=TRUE)
    # Create time stamp for that parameter set; we append a unique random string
    # to ensure the generation of unique names in the case of fast-running models
    now= paste(format(Sys.time(),"%Y%m%d_%H%M%S"),basename(tempfile(pattern="run_")),sep=".")
    # Generate entry in logfile
    write(x=paste(now,"\t",paste(parameters,collapse="\t"),sep=""),
      file=logfile,ncolumns=1, append=TRUE)
    # Save model outputs
    dir_target= paste(outdir_mcs,now,sep="/")
    ok= file.rename(from=outdir_model, to=dir_target)
    if (!ok)
      stop(paste("Cannot save model output in '",outdir_model,"' to '",dir_target,"'.",sep=""))
    ok= dir.create(outdir_model)
    if (!ok)
      stop(paste("Cannot re-create model output directory '",outdir_model,"'.",sep=""))
  } # End of MCS loop
  return(invisible(NULL))
}

################################################################################

#' Evaluation of Monte-Carlo simulation output
#'
#' The function analyzes the output created by \code{\link{mcs_run}}. It
#' computes goodness-of-fit measures for all tested parameter sets.
#'
#' @param outdir_mcs Name of the directory created by \code{\link{mcs_run}}
#'   containing the output of a Monte-Carlo simulation.
#' @param obs_files A \emph{named} vector of file names. Each of the files contains
#'   a time series of observations against which simulation results are compared.
#' @param obs_colsTime A vector of the same length as \code{obs_files}, specifying
#'   the name(s) of the colum(s) holding time information.
#' @param obs_colsValue A vector of the same length as \code{obs_files}, specifying
#'   the name(s) of the colum(s) holding the observation data.
#' @param sim_file \emph{Base name} of the file containing the simulated time
#'   series. The full file path is assembled automatically.
#'
#' @inheritParams modelError_multiDim
#' @inheritParams visEval
#'
#' @return The function returns an object of class 'mcs_eval'. This is a list
#'   with the following components:
#' \item{tbl}{A data frame that lists the tested parameter
#'   sets together with the corresponding values of the objective function(s).
#'   The columns 'err.obs', 'err.len', 'err.fun', and 'err.val' specify the
#'   name of the observation data set, the number of observations,
#'   the name of the objective function, and the value of the objective function, respectively.}
#' \item{par}{A vector holding the names of the parameters as appearing in the column headers of \code{tbl}.}
#'
#' @note Multiple observation files may be supplied in the \code{obs_files}
#'   argument. This makes it possible to analyze the relation between
#'   the model parameters and the simulation error for different sub-sets of
#'   observations. This is often useful if the method of Monte-Carlo simulation
#'   is applied in order to identify a set of optimum model parameters.
#'   For example, in hydrological model calibration, \code{obs_files} could point
#'   to three files: (1) a time series of all observations, (2) a time series of
#'   observations during winter time and (3) another time series of observations
#'   during low-flow periods. This often facilitates the identification of those
#'   model parameters which are sensitive in certain periods of time only (e.g.
#'   during snow cover or low flow).
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export

mcs_eval = function(
  outdir_mcs,
  obs_files,
  obs_colsTime,
  obs_colsValue,
  obs_colsep="\t",
  obs_timeConv= function(x) { as.POSIXct(strptime(x,
    "%Y-%m-%d %H:%M:%S",tz="UTC"),tz="UTC") },
  sim_file,
  sim_colTime,
  sim_colValue,
  sim_colsep="\t",
  sim_timeConv= function(x) { as.POSIXct(strptime(x,
    "%Y-%m-%d %H:%M:%S",tz="UTC"),tz="UTC") },
  obs_nodata,
  gof_function= function(obs,sim) {
    c(bias= mean(sim-obs) ,rmse= sqrt(mean((sim-obs)^2))) 
  }
) {
  # Check args
  if (!dir.exists(outdir_mcs))
    stop(paste("MCS output directory '",outdir_mcs,"' not found."))
  if (is.null(names(obs_files)) || any(names(obs_files) == ""))
    stop("Missing element names in vector of observation files.")
  if (length(obs_colsTime) != length(obs_files))
    stop("Time column must be specified for all observation files.")
  if (length(obs_colsValue) != length(obs_files))
    stop("Value column must be specified for all observation files.")
  testError= gof_function(obs=1:3, sim=2:4)
  if (is.null(names(testError)) || any(names(testError) == ""))
    stop("Missing element names in result vector of G-O-F function.")

  # Read the mcs logfile
  log= read.table(file=paste(outdir_mcs,"mcs.log",sep="/"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
  parNames= names(log)[2:ncol(log)]

  # Read all observation files
  obs= vector(mode="list")
  for (i in 1:length(obs_files)) {
    obs[[i]]= read_timeSeries(file=obs_files[i], colTime=obs_colsTime[i], colValue=obs_colsValue[i],
      colsep=obs_colsep, timeConv=obs_timeConv, nodata=obs_nodata)
    names(obs)[i]= names(obs_files)[i]
    names(obs[[i]])= c("time","obs")
  }

  # Open tempfile
  tf= tempfile()
  ovect=c("timestamp","err.obs","err.len","err.fun","err.val")
  write(x=ovect,file=tf, ncolumns=length(ovect), sep="\t", append=TRUE)

  # Loop through parameter sets
  for (iset in 1:nrow(log)) {
    print(paste("Set ",iset," of ",nrow(log),"...",sep=""))

    # Read simulated time series
    f= paste(outdir_mcs,log$timestamp[iset],sim_file,sep="/")
    if (!file.exists(f))
      stop(paste("Model output file with expected name '",f,"' not found.",sep=""))
    sim= read_timeSeries(file=f, colTime=sim_colTime, colValue=sim_colValue,
      colsep=sim_colsep, timeConv=sim_timeConv, nodata=NA)
    names(sim)= c("time","sim")

    # Compute error for all GOF functions and all observation series
    for (iobs in 1:length(obs)) {
      gofData= merge(x=obs[[iobs]], y=sim, by="time", all=FALSE)
      if (nrow(gofData) == 0)
        stop(paste("No corresponding data in filtered '",obs_files[iobs],"' and '",sim_file,"'.",sep=""))
      gof= gof_function(obs=gofData$obs, sim=gofData$sim)
      for (igof in 1:length(gof)) {
         ovect= c(log$timestamp[iset],names(obs)[iobs],nrow(gofData),names(gof)[igof],gof[igof])
         write(x=ovect, file=tf, ncolumns=length(ovect), sep="\t", append=TRUE)
      }
    }
  }

  # Read tempfile
  res= read.table(file=tf, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  invisible(file.remove(tf))

  # Merge GOF info and parameter values
  res= merge(x=res, y=log, by="timestamp")

  # Return result object with class attribute
  res= list(tbl=res, par= parNames)
  class(res)= "mcs_eval"
  return(res)
}

################################################################################

#' Creates dotty plots from Monte-Carlo simulation output
#'
#' The function takes an output of \code{\link{mcs_eval}} to create
#' plots of the marginal distribution(s) of the objective function(s) with
#' respect to the individual parameters. These plots are sometimes referred to
#' as 'dotty plots'. They provide a means to visually check parameter
#' sensitivity and identifyability. The plots are saved in a file (pdf format).
#'
#' @param x An object of class 'mcs_eval' created by \code{\link{mcs_eval}}.
#' @param annotate Logical switch to enable/disable the display of some statistical
#'   measures in the plots.
#' @param outfile A name for the output file. The extension '.pdf' is appended
#'   automatically.
#' @param overwrite If \code{TRUE}, an existing output file is silently replaced.
#'
#' @return \code{NULL}.
#'
#' @note If the number of objective functions and/or independend observation data sets is
#'   too large the layout may be sub-optimal or plotting may fail at all.
#'   The function is known to work for 3 functions and 3 observation data sets, at least.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export

mcs_plot_dotty= function(
  x,
  annotate=TRUE,
  outfile="dotty.pdf",
  overwrite=FALSE
) {
  # Check args
  if (class(x) != "mcs_eval")
    stop("Argument 'x' must be an object of class 'mcs_eval'.")
  if ((length(annotate) != 1) || (!is.logical(annotate)))
    stop("Argument 'annotate' must be a logical value.")
  if ((length(outfile) != 1) || (!is.character(outfile)))
    stop("Argument 'outfile' must be a character string.")
  if ((length(overwrite) != 1) || (!is.logical(overwrite)))
    stop("Argument 'overwrite' must be a logical value.")
  if (file.exists(outfile) && (!overwrite))
    stop("Output file '",outfile,"' already exists.",sep="")

  # Warn if plotting might fail
  if (length(unique(x$tbl$err.obs)) > 3)
    warning("Number of observation sets > 3. Plotting might fail.")
  if (length(unique(x$tbl$err.fun)) > 3)
    warning("Number of functions > 3. Plotting might fail.")

  # Create dotty plots
  pdf(file=outfile, paper="a4r", height=15/2.54)
  op <- par(mar=c(4, 4, 1.5, 0.5))
  for (p in x$par) {
    layout(matrix(1:(length(unique(x$tbl$err.obs))*length(unique(x$tbl$err.fun))),
      ncol=length(unique(x$tbl$err.obs)), byrow=TRUE))
    for (f in unique(x$tbl$err.fun)) {
      for (o in unique(x$tbl$err.obs)) {
        inds= which((x$tbl$err.obs==o) & (x$tbl$err.fun==f))
        if (!all(is.finite(x$tbl$err.val[inds])))
          warning(paste("Omitting non-finite results (obs.: '",o,"', function: '",f,"').",sep=""))
        yrng= range(x$tbl$err.val[inds], na.rm=TRUE)
        if (!all(is.finite(yrng)))
          stop(paste("Only non-finite results (obs.: '",o,"', function: '",f,"').",sep=""))
        plot(x=range(x$tbl[inds,p]), y=yrng, type="n", xlab=p, ylab=f)
        points(x=x$tbl[inds,p], y=x$tbl$err.val[inds], pch=20)
        abline(h=0, lty=3)
        mtext(side=3, paste("Obs.: ",o,sep=""), cex=0.7)
        if (annotate) {
          legend("right",bty="n",legend=
            c(paste("max",signif(max(x$tbl$err.val[inds],na.rm=TRUE),3),sep="="),
            paste("min",signif(min(x$tbl$err.val[inds],na.rm=TRUE),3),sep="=")),cex=0.7)
        }
      }
    }
  }
  par(op)
  graphics.off()
}

