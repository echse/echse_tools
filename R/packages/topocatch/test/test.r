rm(list=ls())

library("topocatch")

setwd("/home/dkneis/progress/echse/echse_tools/R/packages/topocatch/test")
fileDEM_raw= tempfile()
fileLUSE= tempfile()
fileFCAP= tempfile()
crit_source_area=0.05e06
minlength_reach=100
fill= TRUE
data(dem,landuse,fieldcapacity)
geogrid.writeAscii(dem, file=fileDEM_raw)
geogrid.writeAscii(landuse, file=fileLUSE)
geogrid.writeAscii(fieldcapacity, file=fileFCAP)

outdir="out"

# Set and check output directory
stat= file.info(outdir)
if (is.na(stat$isdir))
  stop(paste("Directory '",outdir,"' does not exist or is not accessible.",sep=""))
if (!stat$isdir)
  stop(paste("Path '",outdir,"' does not point to a directory.",sep=""))
#if (length(list.files(path=outdir)) > 0)
#  stop(paste("Output directory '",outdir,"' is not empty.",sep=""))

fileDEM=        paste(outdir,"dem_filled.asc",sep="/")
fileDIR=        paste(outdir,"dir.asc",sep="/")
fileACC=        paste(outdir,"acc.asc",sep="/")
fileCTI=        paste(outdir,"cti.asc",sep="/")
fileSHP=        paste(outdir,"net.shp",sep="/")
fileCAT=        paste(outdir,"cat.asc",sep="/")
fileHRU=        paste(outdir,"hru.asc",sep="/")
fileAttrCAT=    paste(outdir,"catAttr.txt",sep="/")
fileAttrRCH=    paste(outdir,"rchAttr.txt",sep="/")
fileObjDecl=    paste(outdir,"objDecl.txt",sep="/")
fileObjLink=    paste(outdir,"objLink.txt",sep="/")

prefix_catch= "cat_"
prefix_node= "node_"


################################################################################

# Fill sinks in DEM
if (fill)
  dem.fill(fileIn=fileDEM_raw, fileOut=fileDEM, ndigits=0,
    replace=TRUE, silent=FALSE)

# Analyze DEM
id_outlet= dem.analyze(fileDEM=fileDEM, fileDIR=fileDIR, fileACC=fileACC,
  fileCTI=fileCTI, fileSHP=fileSHP,
  crit_source_area=crit_source_area, x_inBasin=NULL, y_inBasin=NULL,
  id_field="id", class_field="class", minlength_reach=50,
  classname_reach="rch", classname_minireach="minirch", dz_min=0.1,
  replace=TRUE, silent=FALSE
)
print(paste("Computed ID of outlet:",id_outlet))

# Create input for ECHSE-based hydrological modeling
hydroModelData(fileSHP=fileSHP, fileDEM=fileDEM, fileDIR=fileDIR, fileCAT=fileCAT,
  fileAttrCAT=fileAttrCAT, fileAttrRCH=fileAttrRCH, fileObjDecl=fileObjDecl, fileObjLink=fileObjLink,
  id_outlet=id_outlet,
  id_field="id", class_field="class", classname_reach="rch", classname_minireach="minirch",
  classname_node="node", classname_catch="cat", classes_with_catchment= c("rch"),
  nbuffer=1, coord_tol=1.e-09, prefix_node=prefix_node, prefix_catch=prefix_catch,
  updateSHP=TRUE,
  namesIO= data.frame(target=c("qi_avg","qi_end"),source=c("qx_avg","qx_end")),
  replace=TRUE, silent=FALSE
)

# Add further fields to the catchments' attribute table
catAttr= read.table(file=fileAttrCAT, header=TRUE)
catGrid= geogrid.readAscii(fileCAT)

# Add DEM statistics
tmp= geogrid.zones.continuous(grid_zones=catGrid, grid_data=geogrid.readAscii(fileDEM),
  minCoverage=0.95, prefix="elev_", addQuartiles=TRUE, addExtremes=TRUE, sdigits=4, silent=FALSE)
tmp$id= paste(prefix_catch,tmp$id,sep="")
catAttr= merge(x=catAttr, y=tmp, by.x="object", by.y="id", all=TRUE)
if (any(is.na(catAttr)))
  stop("Failed to add elevation info to catchment attributes table.")

# Add land use statistics
luseGrid= geogrid.readAscii(fileLUSE)
tmp= geogrid.zones.classified(grid_zones=catGrid, grid_data=luseGrid, minCoverage=0.95,
  minShare=0.05, keepClasses=numeric(0), prefix="lu", reportArea=FALSE, ndigits=2, silent=FALSE)
tmp$id= paste(prefix_catch,tmp$id,sep="")
catAttr= merge(x=catAttr, y=tmp, by.x="object", by.y="id", all=TRUE)
if (any(is.na(catAttr)))
  stop("Failed to add land use info to catchment attributes table.")

# Convert continuous soil info into classes
fcapGrid= geogrid.readAscii(fileFCAP)
fcapGrid= geogrid.reclass(grid=fcapGrid, lower=c(0,50,100,150,200), upper=c(50,100,150,200,250),
  new=c(1,2,3,4,5), includeLower=TRUE, includeUpper=TRUE)

# Merge info on land use and soil --> HRUs
hruGrid= luseGrid
hruGrid$matrix= hruGrid$matrix + fcapGrid$matrix/10
hruGrid$matrix[fcapGrid$matrix == fcapGrid$nodata_value]= hruGrid$nodata_value
hruGrid$matrix[luseGrid$matrix == luseGrid$nodata_value]= hruGrid$nodata_value
geogrid.writeAscii(hruGrid, file=fileHRU, replace=TRUE)

# Compute areal shares of the HRUs
tmp= geogrid.zones.classified(grid_zones=catGrid, grid_data=hruGrid,
  minCoverage=0.95, minShare=0.05, keepClasses=numeric(0), prefix="hru", reportArea=FALSE, ndigits=2, silent=FALSE)
tmp$id= paste(prefix_catch,tmp$id,sep="")
catAttr= merge(x=catAttr, y=tmp, by.x="object", by.y="id", all=TRUE)
if (any(is.na(catAttr)))
  stop("Failed to add HRU info to catchment attributes table.")

# Rewrite catchments attribute table
write.table(catAttr, file=fileAttrCAT, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

################################################################################
# Create plots

split.screen(c(2,2))

# DEM
screen(1)
geogrid.plot(geogrid.readAscii(fileDEM), title="DEM", randomColors=FALSE, len=6)
# Sub-basins with elevations
screen(2)
showElev= function(tab) { text(round(tab$elev_avg), x=tab$x, y=tab$y, cex=0.6) }
geogrid.plot(catGrid, title="Sub-basin elevations", randomColors=TRUE,  fun=showElev, tab=catAttr)
# HRUs
screen(3)
geogrid.plot(hruGrid, title="HRUs", randomColors=TRUE)

close.screen(all=TRUE)

