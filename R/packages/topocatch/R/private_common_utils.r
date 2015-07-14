
################################################################################
# Function to check arguments

# use len=NULL to skip the length check (for vectors of unknown length and data.frames)
checkArg= function(arg, len=NULL, type="character") {
  caller= as.character(sys.call(-1)[[1]])
  # Get name of arg
  if (missing(arg))
    stop(paste("Argument to be checked is missing.",sep=""))
  name= deparse(substitute(arg))
  # Check args of this function itself
  if (!is.null(len)) {
    if ((length(len) != 1) || (!is.numeric(len)))
      stop("Argument 'len' must be a scalar integer or NULL. Argument check failed.")
  }
  if (length(type) != 1)
    stop("Argument 'type' must be a string representing a data type. Argument check failed.")
  # Perform actual checks
  lenOK= TRUE
  if (!is.null(len)) {
    if (length(arg) != len)
      lenOK= FALSE
  }
  if (type=="character") {
    typeOK= is.character(arg)
  } else if (type=="numeric") {
    typeOK= is.numeric(arg)    
  } else if (type=="integer") {
    typeOK= is.numeric(arg) & (arg == round(arg))
  } else if (type=="logical") {
    typeOK= is.logical(arg)
  } else if (type=="function") {
    typeOK= is.function(arg)
  } else if (type=="data.frame") {
    typeOK= is.data.frame(arg)
  } else {
    stop(paste("Cannot check argument '",name," of function '",caller,"'. Type '",type,"' is unknown.",sep=""))
  }
  if (!all(c(lenOK,typeOK))) {
    if (!is.null(len)) {
      stop(paste("Argument '",name,"' of function '",caller,"' must be an object of type '",type,
        "' with length ",len,".",sep=""))
    } else {
      stop(paste("Argument '",name,"' of function '",caller,"' must be an object of type '",type,
        "'.",sep=""))
    }
  }
  return(invisible(NULL))
}

################################################################################
# Functions to check file and directory names

checkFileIn= function(f) {
  caller= as.character(sys.call(-1)[[1]])
  # Get name of arg
  if (missing(f))
    stop(paste("Argument to be checked is missing.",sep=""))
  name= deparse(substitute(f))
  if ((!is.character(f)) || (length(f) != 1))
    stop(paste("File name in argument '",name,"' of function '",caller,"' must be a character string.",sep=""))
  if (!file.exists(f))
    stop(paste("File '",f,"' supplied as argument '",name,"' of function '",caller,"' does not exist.",sep=""))
}

checkFileOut= function(f, replace=FALSE) {
  caller= as.character(sys.call(-1)[[1]])
  # Get name of arg
  if (missing(f))
    stop(paste("Argument to be checked is missing.",sep=""))
  name= deparse(substitute(f))
  if ((!is.character(f)) || (length(f) != 1))
    stop(paste("File name in argument '",name,"' of function '",caller,"' must be a character string.",sep=""))
  if ((length(replace) != 1) || (!is.logical(replace)))
    replace=FALSE
  if (file.exists(f) && (!replace))
    stop(paste("File '",f,"' supplied as argument '",name,"' of function '",caller,"' already exists.",sep=""))
}

checkDir= function(f, mustBeEmpty=TRUE) {
  caller= as.character(sys.call(-1)[[1]])
  # Get name of arg
  if (missing(f))
    stop(paste("Argument to be checked is missing.",sep=""))
  name= deparse(substitute(f))
  if ((!is.character(f)) || (length(f) != 1))
    stop(paste("Directory name in argument '",name,"' of function '",caller,"' must be a character string.",sep=""))
  stat= file.info(f)
  if (is.na(stat$isdir))
    stop(paste("Directory '",f,"' supplied as argument '",name,"' of function '",caller,"' does not exist or is not accessible.",sep=""))
  if (!stat$isdir)
    stop(paste("Directory '",f,"' supplied as argument '",name,"' of function '",caller,"' is not a valid directory.",sep=""))
  if ((length(mustBeEmpty) != 1) || (!is.logical(mustBeEmpty)))
    mustBeEmpty=FALSE
  if (mustBeEmpty && (length(list.files(path=f)) > 0))
    stop(paste("Directory '",f,"' supplied as argument '",name,"' of function '",caller,"' is not empty.",sep=""))
}


