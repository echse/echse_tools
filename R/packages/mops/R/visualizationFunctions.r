

################################################################################

#' Simple visualization function for manual model calibration
#'
#' The function plots the simulation results for selected data groups and
#' periods of time. Scaling of the plots is determined by the observation data
#' only to facilitate visual comparison with alternative simulation outputs.
#'
#' @param sim_file Name of the file containing the simulated time series.
#' @param sim_colTime Name of column in \code{sim_file} containing time info.
#' @param sim_colValue Name of column in \code{sim_file} containing the data.
#' @param obs_file Name of the file containing the time series of observations.
#' @param obs_colTime Name of column in \code{obs_file} containing time info.
#' @param obs_colValue Name of column in \code{obs_file} containing the data.
#' @param obs_colsep Column separator used in \code{obs_file}. Defaults to TAB.
#' @param obs_timeConv A function to convert the strings in column
#'   \code{obs_colTime} of \code{obs_file}. See argument \code{timeConv} of
#'   the \code{\link{read_timeSeries}} method for details.
#' @param obs_nodata A value which indicates missing or corrupt data in the time
#'   series of observations. Often used values are -9999 or \code{NA}.
#' @param selection_file Name of a file containing the specification of periods
#'   of interest. The file should contain three columns named 'id', 'from', and
#'   'until'. Identifiers for the data selections are expected in the 'id'
#'   field. Those values are treated as character strings. The entries in the
#'   'from' and 'until' columns control which records of the observed and
#'   simulated data contribute to a particular data selection. The entries in 
#'   these fields must be convertable to class POSIXct (see
#'   \code{selection_timeConv}).
#'   A data selection may comprise multiple separate periods of time, thus the
#'   total number of data selections is equal to the number of unique values in
#'   the 'id' field.
#' @param selection_colsep The column separator used in \code{selection_file}.
#' @param selection_timeConv A function to convert the data in the 'from' and
#'   'until' columns of \code{selection_file} into values of class
#'   \code{POSIXct}. See argument \code{obs_timeConv} of method
#'   \code{prefix{modelEvaluationTable}} for details.
#' @param fac_ymin A factor by which the minimum of the observed data is
#'   multiplied when setting the lower limit of the y-axis. Default is 0.
#' @param fac_ymax Like \code{fac_ymax} but for the upper limit of the y-axis.
#'   Default is 1 + 1/4.
#' @param prefix A string used as the initial part of the output file name.
#'
#' @inheritParams modelError_multiDim
#'
#' @return Name of a PDF file for the created plots. The file name is generated
#'   automatically from the system time and the string in \code{prefix}. The
#'   file is created in the current working directory.
#'
#' @note Although one can use this function to plot data of different resolution
#'   (hourly simulation and daily observations, for example), the computed
#'   goodness-of-fit is based on records with identical time stamps (only).
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export

visEval = function(
  sim_file,
  sim_colTime,
  sim_colValue,
  sim_colsep="\t",
  sim_timeConv= function(x) { as.POSIXct(strptime(x,
    "%Y-%m-%d %H:%M:%S",tz="UTC"),tz="UTC") },
  obs_file,
  obs_colTime,
  obs_colValue,
  obs_colsep="\t",
  obs_timeConv= function(x) { as.POSIXct(strptime(x,
    "%Y-%m-%d %H:%M:%S",tz="UTC"),tz="UTC") },
  obs_nodata= -9999,
  selection_file,
  selection_colsep="\t",
  selection_timeConv= function(x) { as.POSIXct(strptime(x,
    "%Y-%m-%d %H:%M:%S",tz="UTC"),tz="UTC") },
  fac_ymin= 0,
  fac_ymax= 1.25,
  prefix= "visEval"
) {
  # Constants
  colNames_sel= list(id= "id", from= "from", until= "until")
  # Check args
  if (!file.exists(obs_file))
    stop(paste("File with observations '",obs_file,"' not found.",sep=""))
  if (!file.exists(sim_file))
    stop(paste("File with simulated data '",sim_file,"' not found.",sep=""))
  if (!file.exists(selection_file))
    stop(paste("File with data selection '",selection_file,"' not found.",sep=""))
  # Read obs data
  obs= read.table(file=obs_file, header=TRUE, sep=obs_colsep)
  if (!all(c(obs_colTime,obs_colValue) %in% names(obs)))
    stop(paste("Missing column(s) in file '",obs_file,"'. Expected colums '",
      obs_colTime,"' and '",obs_colValue,"'.",sep=""))
  obs[,obs_colTime]= obs_timeConv(obs[,obs_colTime])
  obs[(obs[,obs_colValue] == obs_nodata),obs_colValue]= NA
  # Read sim data
  sim= read.table(file=sim_file, header=TRUE, sep=sim_colsep)
  if (!all(c(sim_colTime,sim_colValue) %in% names(sim)))
    stop(paste("Missing column(s) in file '",sim_file,"'. Expected colums '",
      sim_colTime,"' and '",sim_colValue,"'.",sep=""))
  sim[,sim_colTime]= sim_timeConv(sim[,sim_colTime])
  # Open output file
  ofile= paste(prefix,format(Sys.time(),"%Y%m%d%H%M%S"),".pdf",sep="")
  if (file.exists(ofile))
    stop(paste("Output file '",ofile,"' already exists.",sep=""))
  pdf(file=ofile, width=26/2.54, height=18/2.54, onefile=TRUE, paper="a4r")
  # Read file with data selection
  sel= read.table(file=selection_file, header=TRUE, sep=selection_colsep)
  if (!all(colNames_sel %in% names(sel)))
    stop(paste("Missing column names in file with data selection '",selection_file,
    "'. Expecting columns '",paste(colNames_sel,collapse="', '"),"'.",sep=""))
  sel[,colNames_sel$id]= as.character(sel[,colNames_sel$id])
  sel[,colNames_sel$from]= selection_timeConv(sel[,colNames_sel$from])
  sel[,colNames_sel$until]= selection_timeConv(sel[,colNames_sel$until])
  if (any(is.na(sel[,colNames_sel$from])))
    stop(paste("Failed to convert data in column '",colNames_sel$from,"' of file '",
      selection_file,"' to class POSIXct.",sep=""))
  if (any(is.na(sel[,colNames_sel$until])))
    stop(paste("Failed to convert data in column '",colNames_sel$until,"' of file '",
      selection_file,"' to class POSIXct.",sep=""))
  # Process groups
  obs$selected= NA
  sim$selected= NA
  # Loop through groups with unique ID
  for (id in unique(sel[,colNames_sel$id])) {
    obs$selected= FALSE
    sim$selected= FALSE
    print(paste("Processing group '",id,"'...",sep=""))
    sel_thisGroup= sel[sel[,colNames_sel$id] == id,]
    # Select data for all time windows of this group
    for (k in 1:nrow(sel_thisGroup)) {
      if (sel_thisGroup[k,colNames_sel$from] > sel_thisGroup[k,colNames_sel$until])
        stop(paste("Bad time window specified in '",selection_file,
          "'. Error in record ",k," of selection '",id,"'.",sep=""))
      obs$selected[(obs[,obs_colTime] >= sel_thisGroup[k,colNames_sel$from]) & 
        (obs[,obs_colTime] <= sel_thisGroup[k,colNames_sel$until])]= TRUE
      sim$selected[(sim[,sim_colTime] >= sel_thisGroup[k,colNames_sel$from]) & 
        (sim[,sim_colTime] <= sel_thisGroup[k,colNames_sel$until])]= TRUE
    }
    # Create plot
    if (sum(!is.na(obs[obs$selected,obs_colValue])) == 0)
      stop(paste("Found no observed data for selection '",id,"'.",sep=""))
    tlims= range(c(obs[obs$selected,obs_colTime],sim[sim$selected,sim_colTime]))
    ymin= min(obs[obs$selected,obs_colValue], na.rm=TRUE) * fac_ymin
    ymax= max(obs[obs$selected,obs_colValue], na.rm=TRUE) * fac_ymax
    plot(tlims, c(ymin, ymax), type="n",
      xlab=paste("[",min(tlims)," - ",max(tlims),"]"), ylab="Variable")
    lines(obs[obs$selected,obs_colTime], obs[obs$selected,obs_colValue], col="grey")
    lines(sim[sim$selected,sim_colTime], sim[sim$selected,sim_colValue], col="black")
    legend("topright",bty="n",lty=c(1,1),col=c("grey","black"),
      legend=c("obs","sim"))
    legend("bottomright",bty="n",legend=paste("ID: '",id,"'",sep=""))
    # Compute goodness-of-fit based on corresponding datetimes
    gofData= merge(x=obs[obs$selected,c(obs_colTime,obs_colValue)],
      y=sim[sim$selected,c(sim_colTime,sim_colValue)],
      by.x=obs_colTime, by.y=sim_colTime, all=FALSE)
    gofData= gofData[!is.na(gofData[,obs_colValue]),]
    if (nrow(gofData) == 0) {
      txt= paste("Corresp. values: ",nrow(gofData),"    pBias: ","?","%    NS: ",
        "?",sep="")
    } else {
      pbias= (sum(gofData[,sim_colValue])-sum(gofData[,obs_colValue]))/
        sum(gofData[,obs_colValue])*100
      nash= 1-mean((gofData[,sim_colValue] - gofData[,obs_colValue])^2)/
        var(gofData[,obs_colValue])
      txt= paste("Corresp. values: ",nrow(gofData),"    pBias: ",round(pbias,1),
        "%    NS: ",round(nash,2),sep="")
    }
    mtext(side=3,txt)
  }
  graphics.off()
  return(ofile)
}

################################################################################

# Helper function: Statistics for a category
catStat= function(d, cat, gof_function, e.names) {
  tf= tempfile()
  d$level= as.character(format(d$time,format=cat))
  levels= sort(unique(d$level)) 
  for (i in 1:length(levels)) {
    inds= which(d$level == levels[i])
    if (length(inds) > 0) {
      e= gof_function(obs=d$obs[inds], sim=d$sim[inds])
      if (i == 1)
        write(x=paste(names(cat),paste(e.names,collapse="\t"),sep="\t"), file=tf,
          ncolumns=1,append=FALSE)
      write(x=paste(levels[i],paste(e,collapse="\t"),sep="\t"), file=tf,
        ncolumns=1,append=TRUE)
    }
  }
  gof= read.table(file=tf, header=TRUE, sep="\t",
    colClasses= c("character",rep("numeric",length(e))))
  file.remove(tf)
  return(gof)
}

#' Temporal analysis of model error
#'
#' The function analyzes the error of a simulation separately for temporal
#' subsets, currently years and months.
#'
#' @inheritParams modelError_multiDim
#' @inheritParams visEval
#'
#' @return A list of data frames holding the computed results.
#'
#' @note Although one can use this function to plot data of different resolution
#'   (hourly simulation and daily observations, for example), the computed
#'   goodness-of-fit is based on records with identical time stamps (only).
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export

errorDynamics = function(
  obs_file,
  obs_colTime,
  obs_colValue,
  obs_colsep="\t",
  obs_timeConv= function(x) { as.POSIXct(strptime(x,
    "%Y-%m-%d %H:%M:%S",tz="UTC"),tz="UTC") },
  sim_file,
  sim_colTime,
  sim_colValue,
  sim_colsep="\t",
  sim_timeConv= function(x) { as.POSIXct(strptime(x,
    "%Y-%m-%d %H:%M:%S",tz="UTC"),tz="UTC") },
  obs_nodata=-9999,
  periods.select= data.frame(begin=c(), end=c()),
  periods.ignore= data.frame(begin=c(), end=c()),
  gof_function= function(obs,sim) {c(
    nPairs= length(obs),
    pBias= (sum(sim)-sum(obs))/sum(obs)*100,
    NSE= 1-mean((obs-sim)^2)/var(obs)
  )}
) {
  # Check args
  if (!file.exists(obs_file))
    stop(paste("File with observations '",obs_file,"' not found.",sep=""))
  if (!file.exists(sim_file))
    stop(paste("File with simulated data '",sim_file,"' not found.",sep=""))
  # Find corresponding data
  obs= read_timeSeries(file=obs_file, colTime=obs_colTime, colValue=obs_colValue,
    colsep=obs_colsep, timeConv= obs_timeConv, nodata=obs_nodata,
    periods.select=periods.select, periods.ignore=periods.ignore) 
  names(obs)[names(obs) == "value"]= "obs"
  # Read sim data
  sim= read_timeSeries(file=sim_file, colTime=sim_colTime, colValue=sim_colValue,
    colsep=sim_colsep, timeConv= sim_timeConv, nodata=NA,
    periods.select=periods.select, periods.ignore=periods.ignore)
  names(sim)[names(sim) == "value"]= "sim"
  # Merge obs and sim data
  d= merge(x=obs, y=sim, by="time", all=FALSE)
  if (nrow(d) == 0)
    stop(paste("No corresponding data in filtered '",obs_file,"' and '",sim_file,"'.",sep=""))
  rm(obs, sim)
  # Identify/set function names
  e= gof_function(obs=1:10, sim=2:11)  
  if (is.null(names(e))) {
    names(e)= paste("gof",1:length(e),sep="_")
  } else {
    noname= which(names(e) == "")
    if (length(noname) > 0)
      names(e)[noname]= paste("gof",noname,sep="_")
  }
  e.names=names(e)
  rm(e)
  # Set categories
  cats= c(Year="%Y", Month="%m", YearMonth="%Y-%m")
  # Init result
  result= vector(mode="list", length=length(cats))
  # Statistics for categories
  for (i in 1:length(cats)) {
    result[[i]]= catStat(d=d, cat=cats[i], gof_function=gof_function, e.names=e.names)
    names(result)[i]= names(cats)[i]
  }
  return(result)
}


