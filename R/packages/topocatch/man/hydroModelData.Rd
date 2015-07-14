\name{hydroModelData}
\alias{hydroModelData}
\title{Derive input for hydrological catchment modeling from pre-processed geo data}
\usage{
  hydroModelData(fileSHP, fileDEM, fileDIR, fileCAT,
    fileAttrCAT, fileAttrRCH, fileObjDecl, fileObjLink,
    id_outlet, id_field = "id", class_field = "class",
    classname_reach = "rch",
    classname_minireach = "minirch",
    classname_node = "node", classname_catch = "cat",
    classname_gage = "gage",
    classes_with_catchment = c(classname_reach),
    nbuffer = 1, coord_tol = 0.001, prefix_node = "node_",
    prefix_catch = "cat_", min_slope = 1e-04,
    updateSHP = FALSE,
    namesIO = data.frame(target = c("qi"), source = c("qx")),
    gageLocations = data.frame(id = c(), x = c(), y = c()),
    findAffectedGages = FALSE, replace = FALSE,
    silent = TRUE)
}
\arguments{
  \item{fileSHP}{Name/path of an INPUT shape file with line
  features representing drainage lines. This file is either
  an output of \code{\link{dem.analyze}}, a manually edited
  version of such output or an external file (possibly
  created by manual digitizing). If \code{updateSHP} is
  \code{TRUE}, some new fields will be appended to the
  shape file's attribute table (making it an input
  \emph{and} output).}

  \item{fileDEM}{Name/path of the INPUT file containing the
  \emph{sink-filled} DEM (ASCII grid). It is used here to
  estimate the bed slope of river reaches.  If the shape
  file supplied as \code{fileSHP} was created by a call to
  \code{\link{dem.analyze}}, the same DEM as in this call
  should be used.}

  \item{fileDIR}{Name/path of the INPUT file containing
  flow direction codes computed by
  \code{\link{dem.analyze}} (ASCII grid).}

  \item{fileCAT}{Name/path of the OUTPUT grid file showing
  shape, position and extent of the created catchments
  (ASCII grid). Each catchment is identified by an ID
  number. This integer code is derived from the IDs of the
  corresponding features in the input shape file
  (\code{fileSHP}).}

  \item{fileAttrCAT}{Name/path of a tabular OUTPUT file
  listing basic attributes of the generated catchments such
  as the areal extent (in field 'area') and the positions
  of the center of gravity ('x' and 'y' fields). The unit
  of the area is square of the basic length unit used in
  the input files. See also the \code{findAffectedGages}
  argument.}

  \item{fileAttrRCH}{Name/path of a tabular OUTPUT file
  listing basic attributes of reach objects (incl.
  'mini-reaches'). The attributes include the reach length
  (field 'length'), the elevation of the two reach ends
  (fields 'elev_min' and 'elev_max'), a rough estimate of
  the bed slope derived from the former information
  ('slope' field) as well as the total area of the reach's
  upstream catchment (in field 'upstreamArea'). See also
  the \code{findAffectedGages} argument.}

  \item{fileObjDecl}{Name/path of a tabular OUTPUT file to
  be used as an input by ECHSE-based hydrological models.
  This file represents the so-called 'object declaration
  table' holding IDs and class info for all objects.}

  \item{fileObjLink}{Name/path of a tabular OUTPUT file to
  be used as an input by ECHSE-based hydrological models.
  This file represents the so-called 'object linkage
  table'. It provides info on object interactions.}

  \item{id_outlet}{ID of the line feature representing the
  system's outlet (ID of the most downstream reach). The
  value must exist in the ID field of the shape file's
  attribute table specified by \code{id_field}. If the
  shape file supplied as \code{fileSHP} was created by a
  call to \code{\link{dem.analyze}}, without subsequent
  modification, the return value of this function is an
  appropriate input.}

  \item{id_field}{Name of the field in the shape file's
  attribute table containing feature IDs.}

  \item{class_field}{Name of the field in the shape file's
  attribute table holding information on the the features'
  classes. This field is used, for example, to distinguish
  between reaches, pipes, reservoirs, etc. See the
  \code{classes_with_catchment} argument for further
  details.}

  \item{classname_reach}{Class name used for reach objects.
  If the shape file supplied as \code{fileSHP} was created
  by a call to \code{\link{dem.analyze}}, the value should
  be consistent with the one used in that call.}

  \item{classname_minireach}{Class name used for very short
  reach objects with negligible travel time.  If the shape
  file supplied as \code{fileSHP} was created by a call to
  \code{\link{dem.analyze}}, the value should be consistent
  with the one used in that call.}

  \item{classname_node}{Class name to be assigned to node
  objects (junctions).}

  \item{classname_catch}{Class name to be assigned to
  catchment objects.}

  \item{classname_gage}{Class name to be assigned to gage
  objects, if existent.}

  \item{classes_with_catchment}{A vector of class names
  (strings). Catchments are generated only for features
  belonging to those classes. These class names are
  expected to exist in the attribute table's field
  specified by \code{class_field}.  Typically, catchments
  are to be generated for objects of the reach class, at
  least. Lake and reservoir objects, usually have a
  catchment on its own as well.}

  \item{nbuffer}{An integer value to control the thickness
  of the rasterized lines during vector-to-raster
  conversion. The converted lines will be \eqn{nbuffer * 2
  + 1} cells wide. Thus, the default value of 1 results in
  rasterized lines being 3 cells wide.}

  \item{coord_tol}{A small distance value (same spatial
  unit as used in the shape file) to account for precision
  problems. The positions of two lines' end points are
  assumed to be identical if their distance is <= the value
  of \code{coord_tol}. The default is 1.e-03.  Shape files
  created by a call to \code{\link{dem.analyze}} are always
  \emph{exact} in this respect and the small default value
  is appropriate.  If the shape file was created manually
  in a GIS without the use of a snapping mechanism (which
  is \emph{not} recommended), a larger value might be
  required.  The value of \code{coord_tol} is also used
  when gage objects are specified by position via the
  \code{gageLocations} argument.}

  \item{prefix_node}{String used as a prefix when
  generating the IDs of node objects based on the features
  ID of the respective downstream objects.}

  \item{prefix_catch}{String used as a prefix when
  generating the IDs of catchment objects based on the IDs
  of the corresponding features from the shape file.}

  \item{min_slope}{Minimum slope of river bed
  (dimensionless). This value is used as a substitute for
  zero-slope values in the reach attributes table.  Such
  zero-values typically exist in nearly flat areas.}

  \item{updateSHP}{If \code{TRUE}, fields will be appended
  to the shape file's attribute table. Those fields contain
  information for reach-like objects such as the length and
  the total upstream area. The value of -9999 is used in
  these fields for non-reach-like objects and those objects
  which were dropped from the system because of missing
  inflow from upstream.}

  \item{namesIO}{A data frame with 2 columns named 'target'
  and 'source'. Each record represents the names of a pair
  of corresponding input-output variables. A typical
  example would be a data frame with just a single record
  and the strings "inflow" and "outflow" in the 'target'
  and 'source' column, respectively.  Additional records in
  the data frame declare additional I/O variables whose
  values are exchanged between interacting objects.}

  \item{gageLocations}{A data frame with (at least) 3
  columns named 'id', 'x' and 'y'. Each record defines the
  position of a gage. By default, this data frame is empty.
  The positions given by x and y \emph{must} coindide with
  the coordinates of the start or end point of a line
  feature in the shape file. Thus, a gage can only be
  located at either end of a reach but not somewhere in the
  mid. Note that the coordinates must match to the
  precision defined by \code{coord_tol}.}

  \item{findAffectedGages}{If \code{TRUE}, the function
  identifies the gage objects being affected by each
  particular object. For each gage, an additional column is
  appended to the output files \code{fileObjDecl},
  \code{fileAttrCAT}, and \code{fileAttrRCH}. Column names
  are taken from the 'id' field of \code{gageLocations}.
  The values in those columns are 1 or 0, where 1 means
  that the gage (column) is affected by the simulated
  object (row) and 0 indicates the opposite (no
  interaction). The default for \code{findAffectedGages} is
  \code{FALSE} because the current algorithm may be slow
  for very large systems.  If \code{gageLocations} is empty
  this argument setting has no effect.}

  \item{replace}{Should existing output files be silently
  replaced?}

  \item{silent}{Logical value to turn diagnostic messages
  on/off.}
}
\value{
  The function returns \code{NULL}. All computed
  information is written to the respective output files.
}
\description{
  The function identifies the relevant objects for
  object-based hydrological catchment modeling from
  pre-processed geo data sets. It creates several output
  files which can directly serve as an input for
  hydrological models, namely those build with the ECHSE
  simulation environment.
}
\note{
  Not all of the features in the shape file are retained as
  reach-like objects. Those features that do not have a
  catchment upstream are dropped automatically during
  creation of the various output tables. This is true in
  particular for head-reaches.

  If gage objects are to be considered, an iterative
  approach may be convienent.  In the 1st step, this
  function is called without any gage specifications
  (default for \code{gageLocations}). Suitable gage
  positions can then be identified in the GIS based on the
  updated shape file if \code{updateSHP} was set to
  \code{TRUE} (coordinates are available as attributes and
  can be queried). In the 2nd step, the function is called
  with gage positions supplied in \code{gageLocations}.
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}

