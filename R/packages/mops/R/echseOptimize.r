
#' Optimization with ECHSE-based models
#'
#' The function evaluates and optionally minimizes an objective function
#' that compares the output of an ECHSE-based model (single variable
#' at single object) with observation data. The optimized parameters are
#' dimensionless multipliers applied to the model's input.
#'
#' @param opt_func A wrapper for the function that carries out the actual
#'   optimization. It is called as \code{opt_func(f, opt_args)} with \code{f}
#'   being the objective function and \code{opt_args} containing all other
#'   inputs, typically in form of a list. The objective function \code{f} is
#'   expected to take the vector of parameters (i.e. multipliers) as the only argument.
#'   
#' @param opt_args Data to be passed as the second argument to \code{opt_func}.
#'   The necessary contents of this argument depends on the nature of \code{opt_func}.
#' 
#' @param file_cnf Configuration file of the ECHSE-based model. The keywords
#'   'simStart', 'simEnd', and 'outputDirectory' must \emph{NOT} be present in
#'   the file. The same is true for the keywords specified in column 'config_key'
#'   of \code{filesTable} and the \code{model_args} argument.
#'   This is due to the fact that all these values are
#'   automatically passed to the model as command line arguments.
#' @param simStart Start time of the simulation as a value of class \code{POSIXct}.
#' @param simEnd End time of the simulation as a value of class \code{POSICct}.
#' @param outdir Name/path of a non-existing directory where all outputs are to
#'   be stored. An error is generated if the attempt to create this folder fails.
#' @param placeholder_openChar Single character to mark the begin of a
#'   placeholder in the files listed in column 'file_template' of \code{filesTable}.
#' @param placeholder_closeChar Single character to mark the end of a
#'   placeholder in the files listed in column 'file_template' of \code{filesTable}.
#' @param colsep_input Single character used to separate columns in the model's
#'   input files, including the files listed in column 'file_base' of \code{filesTable}.
#' @param filesTable A data frame with columns 'config_key', 'file_template', and 'file_base'.
#'   Entries in column 'config_key' must be valid keywords used in the
#'   configuration data (file or command line) of an ECHSE-based model.
#'   Entries in column 'file_template' represent the names of template files
#'   to be updated via the \code{\link{update_template}} method. These files must
#'   contain placeholders whose names correspond to the optimized multipliers (
#'   passed via \code{opt_args}.
#'   Entries in column 'file_base'  must be the names of model input files which are to
#'   be multiplied by the updated templates (using \code{\link{multiplyFiles}})
#'   before they are actually read by the model. These files must contain
#'   numeric data rather than placeholders but, apart from that, their formats must
#'   be identical to the formats of the corresponding template files.
#' @param sim_objectNames Character string representing the objects of interest.
#'   The model error is analyzed using the model's output for these objects.
#'   In the case of classical optimization, this vector will be of length one
#'   so that \code{func_err} yields a scalar result.
#' @param evalStart A value of class \code{POSIXct} defining the start of the
#'   time period for which the model error is to be analyzed.
#' @param evalEnd A value of class \code{POSIXct} defining the end of the
#'   time period for which the model error is to be analyzed.
#' @param func_err A function with interface \code{func_err(obs, sim)} taking as
#'   input two vectors of observed and simulated values, respectively. It must
#'   return a scalar measure of the model error if \code{opt_func} carries out
#'   function minimization.
#'
#' @inheritParams modelError_multiDim
#' @inheritParams run_model
#' @inheritParams visEval
#'
#' @return The function returns the result value of \code{opt_func}.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export

echse_optimize= function(
  opt_func,
  opt_args,
  model_path,
  model_args,
  file_cnf,
  simStart,
  simEnd,
  outdir,
  placeholder_openChar="{",
  placeholder_closeChar="}",
  colsep_input="\t",
  filesTable,
  sim_objectNames,
  sim_colValue,
  sim_colsep="\t",
  obs_file,
  obs_colTime,
  obs_colsep="\t",
  obs_timeConv= function(x) { as.POSIXct(strptime(x,
    "%Y-%m-%d %H:%M:%S",tz="UTC"),tz="UTC") },
  obs_nodata= -9999,
  evalStart= simStart,
  evalEnd= simEnd,
  func_err= function(obs, sim) {sqrt(mean((obs-sim)^2))},
  unlist=FALSE,
  stopIfNoData=TRUE
) {
  # Check args (further checks done by sub-routines)
  tryCatch({
    stopifnot(is.function(opt_func))
    stopifnot(is.list(opt_args))
    stopifnot(is.character(model_path), length(model_path)==1)
    stopifnot(is.character(file_cnf), length(file_cnf)==1, file.exists(file_cnf))
    stopifnot(class(simStart)[1]=="POSIXct", class(simEnd)[1]=="POSIXct")
    stopifnot(is.character(outdir), length(outdir)==1)
    stopifnot(is.character(placeholder_openChar), length(placeholder_openChar)==1,
      nchar(placeholder_openChar)==1)
    stopifnot(is.character(placeholder_closeChar), length(placeholder_closeChar)==1,
      nchar(placeholder_closeChar)==1)
    stopifnot(is.character(colsep_input), length(colsep_input)==1, nchar(colsep_input)==1)
    stopifnot(is.data.frame(filesTable), nrow(filesTable)>0,
      all(c("config_key","file_template","file_base") %in% names(filesTable)))
    stopifnot(is.character(sim_objectNames), length(sim_objectNames)>=1)
    stopifnot(is.character(sim_colValue), length(sim_colValue)==1)
    stopifnot(is.character(sim_colsep), length(sim_colsep)==1, nchar(sim_colsep)==1)
    stopifnot(is.character(obs_file), length(obs_file)==1)
    stopifnot(is.character(obs_colTime), length(obs_colTime)==1)
    stopifnot(is.character(obs_colsep), length(obs_colsep)==1, nchar(obs_colsep)==1)
    stopifnot(is.function(obs_timeConv), length(obs_timeConv)==1)
    stopifnot(class(evalStart)[1]=="POSIXct")
    stopifnot(class(evalEnd)[1]=="POSIXct")
    stopifnot(is.function(func_err))
    stopifnot(is.logical(unlist), length(unlist)==1)
    stopifnot(is.logical(stopIfNoData), length(stopIfNoData)==1)
  }, error= function(e) {
    stop(paste("Error in arguments. Details: ",e,sep=""))
  })
  # Check output directory
  if (!dir.exists(outdir))
    stop(paste("Output directory '",outdir,"' does not exist.",sep=""))
  if (length(list.files(path=outdir)) != 0)
    stop(paste("Output directory '",outdir,"' is not empty.",sep=""))

  # Read observation data
  obs= read.table(file=obs_file, header=TRUE, sep=obs_colsep)
  obs[,obs_colTime]= obs_timeConv(obs[,obs_colTime])
  if (!is.na(obs_nodata)) {
    for (col in names(obs)[which(names(obs) != obs_colTime)]) {
      obs[(obs[,col] == obs_nodata),col]= NA
    }
  }
  # Choose automatic names for message files
  file_err= paste(outdir,"/",basename(model_path),".err.html",sep="")
  file_log= paste(outdir,"/",basename(model_path),".log",sep="")
  # Generate names for updated files
  filesTable$file_currentFactors= paste(outdir,"/currentFactors__",
    filesTable$config_key,".tmp",sep="")
  filesTable$file_adjustedInputs= paste(outdir,"/adjustedInputs__",
    filesTable$config_key,".tmp",sep="")
  # INTERNAL FUNCTION: Saves files with initial running numbers
  saveOutput= function(p, gof, outdir) {
    file.remove(file_log)
    file.remove(filesTable$file_currentFactors)
    file.remove(filesTable$file_adjustedInputs)

    f= list.files(path=outdir, pattern="^.+[.].+$", full.names=FALSE)
    tmp= substr(f,1,(regexpr(text=f, pattern=".", fixed=TRUE)-1))
    to_be_renamed= (!is.finite(suppressWarnings(as.numeric(tmp))))
    if (sum(to_be_renamed) == length(f)) {
      n=1
    } else {
      n= max(as.numeric(tmp[!to_be_renamed])) + 1
    }
    if (!all(file.rename(from=paste(outdir,"/",f[to_be_renamed],sep=""),
      to=paste(outdir,"/",n,".",f[to_be_renamed],sep=""))))
      stop("Failed to save model output.")

## OLD VERSION (saved files in sub-folders which are difficult to delete from within R)
#    f= list.files(path=outdir, pattern="^.+[.].+$", full.names=FALSE)
#    d= list.dirs(path=outdir, full.names=FALSE, recursive=FALSE)
#    n= length(d)+1
#    if (!dir.create(path=paste(outdir,"/",n,sep=""), showWarnings=TRUE, recursive=FALSE))
#      stop("Failed to create folder to save model output.")
#    if (!all(file.rename(from=paste(outdir,f,sep="/"), to=paste(outdir,n,f,sep="/"))))
#      stop("Failed to save model output.")
## END OLD VERSION

    return(0)
  }
  # INTERNAL: Command line arguments of the model call
  cmdArgs= c(model_args, file_control= file_cnf, file_err= file_err, file_log= file_log,
    format_err= "html", silent= "true",
    simStart= format(simStart,"%Y-%m-%d %H:%M:%S"),
    simEnd= format(simEnd,"%Y-%m-%d %H:%M:%S"),
    outputDirectory= outdir,
    filesTable$file_adjustedInputs)
  names(cmdArgs)[(length(cmdArgs)-nrow(filesTable)+1):length(cmdArgs)]= filesTable$config_key
  # INTERNAL FUNCTION: File update function to be called inside the objective function.
  updateFiles= function(p, argList) {
    stopifnot(all(c("fileTbl","charOpen","charClose","colsep") %in% names(argList)))
    stopifnot(is.data.frame(argList$fileTbl), all(c("file_template","file_base",
      "file_currentFactors","file_adjustedInputs") %in% names(argList$fileTbl)))
   # updating of templates
    for (i in 1:nrow(argList$fileTbl)) {
      n= update_template(file_template=argList$fileTbl$file_template[i],
        vect_updates=p, char_open=argList$charOpen, char_close=argList$charClose,
        file_result=argList$fileTbl$file_currentFactors[i], overwrite=FALSE)
      if (n == 0)
        stop(paste("Updating of '",argList$filesTable$file_template[i],"' failed.",sep=""))
    }
    # multiplication of updated templates with base files
    for (i in 1:nrow(argList$fileTbl)) {
      multiplyFiles(file1=argList$fileTbl$file_base[i],
        file2=argList$fileTbl$file_currentFactors[i],
        outfile=argList$fileTbl$file_adjustedInputs[i], colsep=argList$colsep, overwrite=FALSE)
    }
    return(0)
  }
  # INTERNAL FUNCTION: Wrapper to simplify the objective function
  err= function(parameters) {
    sim_files= paste(outdir,"/",sim_objectNames,".txt",sep="")
    names(sim_files)= sim_objectNames
    return(modelError_multiDim(
      parameters=parameters,
      model_path= model_path,
      model_args= cmdArgs,
      func_first=updateFiles,
      moreArgs_first=list(fileTbl=filesTable, charOpen=placeholder_openChar, 
        charClose=placeholder_closeChar, colsep=colsep_input),
      func_final=saveOutput,
      moreArgs_final=outdir,
      observed= obs,
      obs_colTime= obs_colTime,
      sim_files= sim_files,
      sim_colsTime=rep("end_of_interval", length(sim_files)),
      sim_colsValue=rep(sim_colValue, length(sim_files)),
      sim_colsep=sim_colsep,
      periods.select= data.frame(begin=evalStart, end=evalEnd),
      periods.ignore= data.frame(begin=c(), end=c()),
      gof_function=func_err,
      unlist=unlist,
      stopIfNoData=stopIfNoData
    ))
  }
  # Call optimization routine
  opt= opt_func(f=err, opt_args)
  # Return
  return(opt)
}


