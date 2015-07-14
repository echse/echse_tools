rm(list=ls())
options(warn=2)

# Script for interactive plotting of multiple time series in R.
# david.kneis@uni-potsdam.de, 2008-2012

# Zooming and panning is based on the zoom() function by Huidae Cho (see below).

################################################################################
### PART 1 - FUNCTIONS #########################################################
################################################################################

push_to_clipboard= function(x) {
  if ((!identical(class(x),"character")) || (length(x) > 1))
    stop("Expecting a character string as input.")
  if (.Platform$OS.type == "windows") {
    write(x=x, file="clipboard") 
  } else if (.Platform$OS.type == "unix") {
    res= system(paste("echo -n '",x,"' | xclip -selection clipboard",sep=""))
    if (res != 0) {
      stop(paste("Failed to push string '",x,"' to clipboard. Is 'xclip' installed?",sep=""))
    }
  } else {
    errstop("Platform not supported.")
  }
}

safeLocator= function() {
  x= try(locator(1), silent=TRUE)
  if (class(x) == "try-error") {
    q()
  } else {
    return(x)
  }
}

inBox= function(x, y, xmin, xmax, ymin, ymax) {
  return((x > xmin) & (x < xmax) & (y > ymin) & (y < ymax))
}

zoom = function(plotfunc, extent){
	# Original function: (C) 2007 GPL by Huidae Cho <http://geni.ath.cx>
	lim = extent
	lim.stack = c(lim$x, lim$y)
	lim.depth = 1
	lim.cur = 1
	repeat{

    # Define buttons
    buttons= data.frame(stringsAsFactors=FALSE,
                                label="Full size",colBG="darkgrey",colFG="white")
    buttons= rbind(buttons,list(label="Zoom out X",colBG="grey",colFG="black"))
    buttons= rbind(buttons,list(label="Zoom out Y",colBG="grey",colFG="black"))
    buttons= rbind(buttons,list(label="Zoom out XY",colBG="grey",colFG="black"))
    buttons= rbind(buttons,list(label="Shift Right",colBG="darkgrey",colFG="white"))
    buttons= rbind(buttons,list(label="Shift Left",colBG="darkgrey",colFG="white"))
    buttons= rbind(buttons,list(label="Shift Up",colBG="darkgrey",colFG="white"))
    buttons= rbind(buttons,list(label="Shift Down",colBG="darkgrey",colFG="white"))
    buttons= rbind(buttons,list(label="Previous zoom",colBG="grey",colFG="black"))
    buttons= rbind(buttons,list(label="Next zoom",colBG="grey",colFG="black"))
    buttons= rbind(buttons,list(label="Quit",colBG="darkgrey",colFG="white"))
    buttons$xmin=NA
    buttons$xmax=NA
    buttons$ymin=NA
    buttons$ymax=NA

    # Split the sreen into a matrix
    # Bottons appear in column 1
    # Legend appears in row 1, column 2...nc
    # Plot appears in column 2...nc, row 2...nr
    NROWS=nrow(buttons)+1
    NCOLS=10
    layout(matrix(c(rep(1,NROWS),rep(c(2,rep(3,(NROWS-1))),(NCOLS-1))),ncol=NCOLS,nrow=NROWS,byrow=FALSE))
    par(mar=c(0,0,0,0)+0.05)


    # Draw buttons in column 1 and save corresponding coordinates (NDC)
    plot(0:1,0:1,type="n",xlab="",ylab="",xaxt="n",yaxt="n",xaxs="i",yaxs="i",bty="n")
    coords.fig= par("fig")
    coords.plt= par("plt")
    for (i in 1:nrow(buttons)) {
      # Draw
      xL=0
      xR=1
      yB=1 - 1/nrow(buttons) * i
      yT=1 - 1/nrow(buttons) * (i-1)
      rect(xleft=xL,xright=xR,ybottom=yB,ytop=yT,col=buttons$colBG[i],border="white")
      text(x=mean(c(xL,xR)), y=mean(c(yB,yT)), buttons$label[i], col=buttons$colFG[i])
      # Save coordinates
      buttons$xmin[i]=coords.fig[1]+coords.plt[1]*(coords.fig[2]-coords.fig[1])
      buttons$xmax[i]=coords.fig[1]+coords.plt[2]*(coords.fig[2]-coords.fig[1])
      buttons$ymin[i]=coords.fig[3]+(coords.plt[3]+(nrow(buttons)-i)/nrow(buttons))*(coords.fig[4] - coords.fig[3])
      buttons$ymax[i]=coords.fig[3]+(coords.plt[3]+(nrow(buttons)-i+1)/nrow(buttons))*(coords.fig[4] - coords.fig[3])
    }

    # Draw legend and series
    plot.legend()
		plot.series(lim)
    coords.fig= par("fig")
    coords.plt= par("plt")

    # Zoom loop
		l = safeLocator()
		if(is.null(l))
			break
    abline(v=l$x, lty=3)  # Show first location
    abline(h=l$y, lty=3)  # Show first location
		ext = par()$usr
    # Check for activated buttons
    frac.x= (l$x-ext[1])/(ext[2]-ext[1])
    frac.y= (l$y-ext[3])/(ext[4]-ext[3])   
    x_device= coords.fig[1] + (coords.plt[1]+frac.x*(coords.plt[2]-coords.plt[1])) * (coords.fig[2]-coords.fig[1])
    y_device= coords.fig[3] + (coords.plt[3]+frac.y*(coords.plt[4]-coords.plt[3])) * (coords.fig[4]-coords.fig[3])
    indexButton= which(inBox(
      x=x_device, y=y_device,
      xmin=buttons$xmin, xmax=buttons$xmax,
      ymin=buttons$ymin, ymax=buttons$ymax
    ))

# Use for debug purposes
#    text(x=mean(c(ext[1],ext[2])),y=mean(c(ext[3],ext[4])), paste(
#      "x",x_device,
#      "y",y_device,"\n",
#       "fig.x1", coords.fig[1], "fig.x2",coords.fig[2],
#       "fig.y1", coords.fig[3], "fig.y2",coords.fig[4], "\n",
#       "plt.x1", coords.plt[1], "plt.x2",coords.plt[2],
#       "plt.y1", coords.plt[3], "plt.y2",coords.plt[4]
#    ) )
#    text(x=mean(c(ext[1],ext[2])),y=mean(c(ext[3],ext[4])), paste(
#       buttons$xmin[1],buttons$xmax[1],buttons$ymin[1],buttons$ymax[1],"\n",
#       buttons$xmin[7],buttons$xmax[7],buttons$ymin[7],buttons$ymax[7]
#    ) )


    # No botton - wait for 2nd cursor action (zoom in or re-center)
    if (length(indexButton) == 0) {
		  l2 = safeLocator()
		  if (is.null(l2))
			  break
      # Zoom in
		  if (sum(l$x == l2$x) || sum(l$y == l2$y)) {
			  w = lim$x[2] - lim$x[1]
			  h = lim$y[2] - lim$y[1]
			  lim = list(x=c(l2$x-w/2, l2$x+w/2),
				  y=c(l2$y-h/2, l2$y+h/2))
      # Center
		  } else {
			  lim = list(x=sort(c(l$x, l2$x)), y=sort(c(l$y, l2$y)))
      }
    # Error
    } else if (length(indexButton) > 1) {
      errstop(paste("Multiple buttons activated.",sep=""))
    # Button action
    } else {
      if (buttons$label[indexButton] == "Quit") {
        break
      } else if ( buttons$label[indexButton] == "Full size") {
  			lim = extent
      } else if ( buttons$label[indexButton] == "Zoom out XY") {
        dx= 0.25*(lim$x[2] - lim$x[1])
        dy= 0.25*(lim$y[2] - lim$y[1])
			  lim = list(x=c(lim$x[1]-dx, lim$x[2]+dx), y=c(lim$y[1]-dy, lim$y[2]+dy))
      } else if ( buttons$label[indexButton] == "Zoom out X") {
        dx= 0.25*(lim$x[2] - lim$x[1])
			  lim = list(x=c(lim$x[1]-dx, lim$x[2]+dx), y=c(lim$y[1], lim$y[2]))
      } else if ( buttons$label[indexButton] == "Zoom out Y") {
        dy= 0.25*(lim$y[2] - lim$y[1])
			  lim = list(x=c(lim$x[1], lim$x[2]), y=c(lim$y[1]-dy, lim$y[2]+dy))
      } else if ( buttons$label[indexButton] == "Shift Right") {
        dx= 0.25*(lim$x[2] - lim$x[1])
			  lim = list(x=c(lim$x[1]+dx, lim$x[2]+dx), y=c(lim$y[1], lim$y[2]))
      } else if ( buttons$label[indexButton] == "Shift Left") {
        dx= 0.25*(lim$x[2] - lim$x[1])
			  lim = list(x=c(lim$x[1]-dx, lim$x[2]-dx), y=c(lim$y[1], lim$y[2]))
      } else if ( buttons$label[indexButton] == "Shift Up") {
        dy= 0.25*(lim$y[2] - lim$y[1])
			  lim = list(x=c(lim$x[1], lim$x[2]), y=c(lim$y[1]+dy, lim$y[2]+dy))
      } else if ( buttons$label[indexButton] == "Shift Down") {
        dy= 0.25*(lim$y[2] - lim$y[1])
			  lim = list(x=c(lim$x[1], lim$x[2]), y=c(lim$y[1]-dy, lim$y[2]-dy))
      } else if ( buttons$label[indexButton] == "Previous zoom") {
	  		cur = lim.cur
        lim.cur= max(lim.cur-1, 1) 
			  if (lim.cur != cur)
				  lim = list(x=lim.stack[lim.cur, 1:2], y=lim.stack[lim.cur, 3:4])
			  next
      } else if ( buttons$label[indexButton] == "Next zoom") {
		  	cur = lim.cur
        lim.cur= min(lim.cur+1, lim.depth)
			  if (lim.cur != cur)
				  lim = list(x=lim.stack[lim.cur, 1:2], y=lim.stack[lim.cur, 3:4])
			  next
      } else {
        errstop(paste("Undefined button activated.",sep=""))
      }
    }
    # Update stack
		if (lim.cur < lim.depth) {
			lim.stack = lim.stack[-((lim.cur+1):lim.depth),]
			lim.depth = lim.cur
		}
		lim.stack = rbind(lim.stack, c(lim$x, lim$y))
		lim.depth = lim.depth + 1
		lim.cur = lim.cur + 1
	} # end of repeat loop
	return(lim)
}

################################################################################

# A function to determine a reasonable time format
timeFmt= function(tmin, tmax) {
  secrange= as.numeric(difftime(tmax, tmin, units="secs"))
  if (secrange > 86400*62) {
    fmt="%d %b '%y"
  } else if ((secrange > 86400*7) && (secrange <= 86400*62)) {
    fmt="%d %b"
  } else if ((secrange > 86400) && (secrange <= 86400*7)) {
    fmt="%d %b %Hh"
  } else {
    fmt="%H:%M:%S"
  }
  return(fmt)
}

################################################################################

# A function to create a vector of times for labelling a time axis. This
# algorithm gives better results than axis.POSIXt.
prettyTimes= function(tmin, tmax, n) {
  # Correct unreasonable n
  if (n < 2) n= 2
  # Years
  s= seq.POSIXt(from= tmin, to=tmax, by="year")
  k= length(s) 
  if (k >= n) {
    i= 1:n * ceiling(k/n)
    t1= ISOdatetime(format(s[1],"%Y"),1,1,0,0,0)
    t2= ISOdatetime(format(s[k],"%Y"),1,1,0,0,0)
    x= seq.POSIXt(from=t1, to=t2, by="year")
    return(list(values=x[i[i<=k]], format=timeFmt(tmin,tmax)))
  }
  # Months
  s= seq.POSIXt(from= tmin, to=tmax, by="month")
  k= length(s) 
  if (k >= n) {
    i= 1:n * ceiling(k/n)
    t1= ISOdatetime(format(s[1],"%Y"),format(s[1],"%m"),1,0,0,0)
    t2= ISOdatetime(format(s[k],"%Y"),format(s[k],"%m"),1,0,0,0)
    x= seq.POSIXt(from=t1, to=t2, by="month")
    return(list(values=x[i[i<=k]], format=timeFmt(tmin,tmax)))
  }
  # Days
  s= seq.POSIXt(from= tmin, to=tmax, by="day")
  k= length(s) 
  if (k >= n) {
    i= 1:n * ceiling(k/n)
    t1= ISOdatetime(format(s[1],"%Y"),format(s[1],"%m"),format(s[1],"%d"),0,0,0)
    t2= ISOdatetime(format(s[k],"%Y"),format(s[k],"%m"),format(s[k],"%d"),0,0,0)
    x= seq.POSIXt(from=t1, to=t2, by="day")
    return(list(values=x[i[i<=k]], format=timeFmt(tmin,tmax)))
  }
  # Hours
  s= seq.POSIXt(from= tmin, to=tmax, by="hour")
  k= length(s) 
  if (k >= n) {
    i= 1:n * ceiling(k/n)
    t1= ISOdatetime(format(s[1],"%Y"),format(s[1],"%m"),format(s[1],"%d"),format(s[1],"%H"),0,0)
    t2= ISOdatetime(format(s[k],"%Y"),format(s[k],"%m"),format(s[k],"%d"),format(s[k],"%H"),0,0)
    x= seq.POSIXt(from=t1, to=t2, by="hour")
    return(list(values=x[i[i<=k]], format=timeFmt(tmin,tmax)))
  }
  # Minutes
  s= seq.POSIXt(from= tmin, to=tmax, by="min")
  k= length(s) 
  if (k >= n) {
    i= 1:n * ceiling(k/n)
    t1= ISOdatetime(format(s[1],"%Y"),format(s[1],"%m"),format(s[1],"%d"),format(s[1],"%H"),format(s[1],"%M"),0)
    t2= ISOdatetime(format(s[k],"%Y"),format(s[k],"%m"),format(s[k],"%d"),format(s[k],"%H"),format(s[k],"%M"),0)
    x= seq.POSIXt(from=t1, to=t2, by="min")
    return(list(values=x[i[i<=k]], format=timeFmt(tmin,tmax)))
  }
  # Seconds
  x= seq.POSIXt(from= tmin, to=tmax, length.out=n)
  return(list(values=x, format=timeFmt(tmin,tmax)))
}

################################################################################


# Function to plot the legend
plot.legend = function() {
par(cex=1)
par(mar=c(0.1,0.1,1,0.1)+0.1)
plot(0:1,0:1,type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
lsym= as.integer(inst$type %in% c("l","S","s"))
psym= (as.integer(inst$type=="p")*20)
psym[psym==0]= NA
legend("center",horiz=TRUE,
  lty=lsym[!is.na(inst$name)],
  pch=psym[!is.na(inst$name)],
  col=inst$color[!is.na(inst$name)],
  legend=inst$name[!is.na(inst$name)],
  bty="n")
}

# Function to plot the time series
plot.series = function(lim) {
  par(cex=1)
  par(mar=c(6,4,1,1)+0.1)
  plot(lim$x,lim$y,type="n",xaxt="n",xlab="",ylab="",xlim=lim$x,ylim=lim$y,xaxs="i",main="")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="lightgrey")
  for (i in 1:nrow(inst)) {
    if (inst$type[i] %in% c("l","S","s")) {
      lines(tser[[i]]$unixtime,tser[[i]]$value,type=inst$type[i],lty=1,
        col=inst$color[i])
    } else if (inst$type[i] == "b") {
      lines(tser[[i]]$unixtime,tser[[i]]$value,lty=1,col=inst$color[i])
      points(tser[[i]]$unixtime,tser[[i]]$value,pch=20,col=inst$color[i])
    } else if (inst$type[i] == "p") {
      points(tser[[i]]$unixtime,tser[[i]]$value,type=inst$type[i],pch=20,
        col=inst$color[i])
    } else {
      errstop(paste("Plot style '",inst$type[i],"' not supported.",sep=""))
    }
  }
  # Set up the time axis (Important: Env. var. TZ must be set properly!)
  tmin= ISOdatetime(1970,1,1,0,0,0,tz="UTC") + as.difftime(lim$x[1], units="secs")
  tmax= ISOdatetime(1970,1,1,0,0,0,tz="UTC") + as.difftime(lim$x[2], units="secs")
  times= prettyTimes(tmin, tmax, 10)
  axis.POSIXct(1, at=times$values, format=times$format, las=2, xlab="")
  mtext(side=3,text=paste(c(tmin,tmax), collapse=" - "), padj=2.5, col="white")
  # Copy current window extent to clipboard
  xleft=  ISOdatetime(1970,1,1,0,0,0,tz="UTC") + as.difftime(par("usr")[1], units="secs")
  xright= ISOdatetime(1970,1,1,0,0,0,tz="UTC") + as.difftime(par("usr")[2], units="secs")
  #ybottom= par("usr")[3]
  #ytop= par("usr")[4]
  #push_to_clipboard(paste(xleft,xright,ybottom,ytop,sep="\t"))
  push_to_clipboard(paste(xleft,xright,sep="\t"))
}

################################################################################

# Opens a new plot window
plotWindow= function(width, height) {
  if (.Platform$OS.type == "windows") {
    windows(width=width,height=height)
  } else if (.Platform$OS.type == "unix") {
    X11(width=width,height=height)
  } else {
    errstop("Platform not supported.")
  }
}

################################################################################

# Read instructions file
readInstructions= function(ifile) {
  if (!file.exists(ifile))
    errstop(paste("Instructions file '",ifile,"' not found.",sep=""))
  inst= read.table(file=ifile,sep="",header=TRUE,comment.char=";",
    blank.lines.skip=TRUE,colClasses=c("character"))
  # Check items in instructions file
  if (nrow(inst) < 1) errstop("Empty instructions file.")
  cols= c("file","colname","color","type","name","header","comment","nodata","factor","xcut")
  if (any(is.na(match(cols,names(inst))))) {
    errstop(paste("Missing column(s) in instructions file.\n",
      "Expected: ",paste(cols,sep="",collapse=", "),".",sep=""))
  }
  inst$nodata= as.numeric(inst$nodata)
  return(inst)
}

################################################################################

# Read all time series and set initial extent (Important: Env. var. TZ must be set properly)
readSeries= function(inst) {
  tser= vector("list",nrow(inst))
  xlims_fixed= FALSE
  xmin= 1e+30
  xmax= -1e+30
  ymin= 1e+30
  ymax= -1e+30
  for (i in 1:nrow(inst)) {
    if (!file.exists(inst$file[i]))
      errstop(paste("Error in record '",i,"' of plot instructions file.",
        " Time series file '",basename(inst$file[i]),"' does not exist.",sep=""))
    tser[[i]]= read.table(file=inst$file[i],sep="\t",header=as.logical(inst$header[i]),
      comment.char=inst$comment[i],blank.lines.skip=TRUE,
      colClasses=c("character"))
    colindex= match(make.names(inst$colname[i]),names(tser[[i]]))
    if (is.na(colindex))
      errstop(paste("Error in record '",i,"' of plot instructions file.",
        "Time series file '",basename(inst$file[i]),"' has no column '",inst$colname[i],"'.",sep=""))
    tser[[i]]= tser[[i]][,c(1,colindex)]
    names(tser[[i]])= c("unixtime","value")
    if (nchar(tser[[i]]$unixtime[1]) == "19") {
      tser[[i]]$unixtime= as.POSIXct(strptime(tser[[i]]$unixtime,"%Y-%m-%d %H:%M:%S",tz="UTC"))
    } else if (nchar(tser[[i]]$unixtime[1]) == "10") {
      tser[[i]]$unixtime= as.POSIXct(strptime(tser[[i]]$unixtime,"%Y-%m-%d",tz="UTC"))
    } else {
      tser[[i]]$unixtime= NA
    }
    tser[[i]]$value= as.numeric(tser[[i]]$value) * as.numeric(inst$factor[i])

    tser[[i]]$unixtime= as.numeric(tser[[i]]$unixtime)
    if (any(is.na(tser[[i]]$unixtime)))
      errstop(paste("Error in array of times in time series file '",
        basename(inst$file[i]),"' specified in record '",i,"' of plot instructions file.",sep=""))
    if (!is.na(inst$nodata[i])) {
      tser[[i]]$value[tser[[i]]$value == inst$nodata[i]]= NA
    } else {
      tser[[i]]$value[is.na(tser[[i]]$value)]= NA
    }    
    # Set x limits
    if (as.logical(inst$xcut[i])) {
      xmin= min(tser[[i]]$unixtime)
      xmax= max(tser[[i]]$unixtime)
      xlims_fixed= TRUE
    } else {
      if (! xlims_fixed) {
        if (min(tser[[i]]$unixtime) < xmin) xmin= min(tser[[i]]$unixtime)
        if (max(tser[[i]]$unixtime) > xmax) xmax= max(tser[[i]]$unixtime)
      }
    }
    # Set initial y limits
    if (min(tser[[i]]$value[!is.na(tser[[i]]$value)]) < ymin) ymin= min(tser[[i]]$value[!is.na(tser[[i]]$value)])
    if (max(tser[[i]]$value[!is.na(tser[[i]]$value)]) > ymax) ymax= max(tser[[i]]$value[!is.na(tser[[i]]$value)])
  }
  return(list(series=tser, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax))
}

################################################################################

# Function to display an error message
breakText= function(text, maxLen) {
  n= 0
  for (i in 1:nchar(text)) {
    n= n + 1
    if (substr(text,i,i) %in% c(" ","\t")) pos=i
    if ((n >= maxLen) && (pos > 0)) {
      text= paste(substr(text,1,(pos-1)),substr(text,(pos+1),nchar(text)),sep="\n")
      n= i-pos
      pos= 0
    }
  }
  return(text)
}

errstop = function (errmsg) {
  plotWindow(width=6,height=3)
  plot(0:1,0:1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
  mtext("Plotting failed. See error message below:",side=3,col="red")
  text(0.5, 0.5, breakText(errmsg,60), cex=0.8)
  mtext("Click on this window to exit.",side=1,col="red")
  l= locator(1)
  q()
}

################################################################################
### PART 2 - PROGRAM ###########################################################
################################################################################

# Get command lines
cmdargs= commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0)
  errstop("Expecting name of instructions file as command line argument.")
ifile= cmdargs[1]

# Make sure that the time zone is set to UTC.
# This is essential for proper READING and PLOTTING of data in its original
# time zone, i.e. to prevent implicit conversion to the current local time zone
# by methods like 'read.table', 'plot', and 'axis.POSIXct'.
Sys.setenv(TZ="UTC")
tz= Sys.getenv("TZ")
if (tz != "UTC") {
  errstop(paste("Environment variable 'TZ' == '",tz,"' but should be 'UTC'."))
}

# Read the instructions file
inst= try(readInstructions(ifile), silent=TRUE)
if (class(inst) == "try-error") errstop(inst)

# Read data and set initial extent
data= try(readSeries(inst), silent=TRUE)
if (class(data) == "try-error") errstop(data)

# Copy series to global variable
tser= data$series

# Set full plot extent
extent= list(x=c(data$xmin,data$xmax),y=c(data$ymin,data$ymax))

# Plot with zoom option
plotWindow(width=12, height=7)
x= try(zoom(plotfunc, extent), silent=TRUE)
if (class(x) == "try-error") errstop(x)

