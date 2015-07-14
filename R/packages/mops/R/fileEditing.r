#' Update a template file containing placeholders
#'
#' The function substitutes placeholders in a template file by actual value(s)
#' provided in a named vector. It can be used in the automatic calibration
#' of models which read numerical parameters from text files.
#'
#' @param file_template Name of the file containing the placeholders to be
#'    updated.
#' @param vect_updates A \emph{named} vector containing the values by which the
#'    placeholders in \code{file_template} are to be substituted. Typically,
#'    this is a numeric vector or a vector of character strings.
#' @param char_open Single character to identify the start of a placeholder in
#'   \code{file_template}.
#' @param char_close Single character to identify the end of a placeholder in
#'   \code{file_template}.
#' @param file_result Name of the output file. The contents is the same as in
#'   \code{file_template} but placeholders are substituted by the matching
#'   values in \code{vect_update}. See details.
#' @param overwrite If \code{FALSE} (default), existing results files will not
#'   be overwritten.
#'
#' @return The function returns an integer value being the number of
#'   placeholders which have been substituted in \code{file_template}.
#'
#' @note For successful updating of the template file, the string between the
#'   opening and closing character of a placeholder must match with the name of
#'   an element in \code{vect_updates}. For example, given that
#'   \code{char_open='<'} and \code{char_open='>'}, a placeholder \emph{<p1>} in
#'   the template file will be substituted by the element of \code{vect_updates}
#'   whose name is 'p1'.
#'
#'   The substitution values are printed using R's default formatting rules. If
#'   a specific number format is required, \code{vect_updates} should be a
#'   vector of character strings holding the properly formatted numbers rather
#'   than a numeric vector.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#'
#' @examples
#' # Create simple template file with placeholders for two values
#' fin= tempfile()
#' write(file=fin, x='p1 {p1} p2 {p2}')
#' # Update the template
#' fout= tempfile()
#' n= update_template(file_template= fin, vect_updates=c(p1=5.0, p2=3.0),
#'   file_result=fout)
#' print(paste("Updated locations:",n))
#' # Show result and clean up
#' x= readLines(con=fout)
#' print(x)
#' file.remove(fin)
#' file.remove(fout)

update_template = function (
  file_template,
  vect_updates,
  char_open="{",
  char_close="}",
  file_result,
  overwrite=FALSE
) {
  # Check arguments
  if (!file.exists(file_template))
    stop(paste("Cannot open template file '",file_template,"'.",sep=""))
  if (length(vect_updates)==0)
    stop("Vector of substitute values has zero length.")
  if (is.null(names(vect_updates)))
    stop("Vector of substitute values is lacking element names.")
  if (any(names(vect_updates) == ""))
    stop("Detected unnamed element(s) in vector of substitute values.")
  if (((nchar(char_open) + nchar(char_close)) != 2) | (char_open == char_close))
    stop("Bad delimiter(s) for placeholders specified.")
  if (file.exists(file_result) && (!overwrite))
    stop(paste("Result file '",file_result,"' already exists.",sep=""))
  # Read the template file
  rec= readLines(file_template)
  numRepl= 0
  # Loop through the lines
  for (i in 1:length(rec)) {
    # Process the current line
    subst= TRUE
    while (subst) {  # There may be multiple placeholders at a single line...
      pos_beg= regexpr(char_open,rec[i],fixed=TRUE)
      pos_end= regexpr(char_close,rec[i],fixed=TRUE)
      if ((pos_beg == -1) && (pos_end == -1)) {  # Nothing to substitute
        subst= FALSE
      } else if ((min(pos_beg,pos_end) == -1) && (max(pos_beg,pos_end) != -1)) {
        stop(paste("Incomplete placeholder at line ",i," of '",
          file_template,"'.",sep=""))
      } else if (pos_beg >= pos_end) { # Bad marker position
        stop(paste("Bad placeholder at line ",i," of '",file_template,"'.",sep=""))
      } else {                         # Substitute substring in markers
        subst= TRUE
        # Get placeholder string
        markstr= substr(rec[i],(pos_beg+1),(pos_end-1))
        if (markstr == "")
          stop(paste("Empty placeholder at line ",i," of '",file_template,"'.",sep=""))
        # Get remaining parts of the current line
        if (pos_beg == 1) {
          first= ""
        } else {
          first= substr(rec[i],1,(pos_beg-1))
        }
        if (pos_end == nchar(rec[i])) {
          final= ""
        } else {
          final= substr(rec[i],(pos_end+1),nchar(rec[i]))
        }
        # Find substitute for the placeholder string in named vector of parameters
        n= which(names(vect_updates) %in% markstr)
        if (length(n) == 0)
          stop(paste("Undefined placeholder in '",char_open,markstr,char_close,
            "' at line ",i," of '",file_template,"'. Name is not present in argument 'vect_updates'.",sep=""))
        if (length(n) > 1)
          stop(paste("Ambiguous placeholder in '",char_open,markstr,char_close,
            "' at line ",i," of '",file_template,"'. Name points to multiple elements of argument 'vect_updates'.",sep=""))
        # Substitute the placeholder string
        rec[i]= paste(first,vect_updates[[n]],final,sep="")
        numRepl= numRepl + 1
      }
    }
  }
  # Write all updated records to output
  for (i in 1:length(rec))
    write(rec[i],file=file_result,append=(i!=1))
  # Return the number of replacements done
  return(numRepl)
}

################################################################################

#' Multiplication of text files of identical layout
#'
#' The function multiplies corresponding values of the \emph{numeric} columns
#' of two text files. The result, including non-numeric columns, is written to
#' a result file.
#'
#' @param file1 First input file containing numeric data in one or more columns.
#'   Column headers are mandatory.
#' @param file2 Second input file with identical layout as \code{file1}. In
#'   particular, column headers, the number of records, and the contents of any
#'   non-numeric columns must be identical.
#' @param outfile Name/path of the result file.
#' @param colsep Character to be used as a column separator in input and output files.
#' @param overwrite Is it OK to replace an already existing output file?

#' @return \code{NULL}.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#'
#' @examples
#' f= tempfile(c("in1","in2","out"))
#' t1= data.frame(city=c("Rome","Berlin","Oslo"),temperature=c(20,10,1))
#' t2= data.frame(city=c("Rome","Berlin","Oslo"),temperature=c(0.5,1,10))
#' write.table(t1, file=f[1], sep="\t", row.names=FALSE, quote=FALSE)
#' write.table(t2, file=f[2], sep="\t", row.names=FALSE, quote=FALSE)
#' multiplyFiles(file1=f[1], file2=f[2], outfile=f[3])
#' t3= read.table(file=f[3], header=TRUE)
#' print(t3)
#' invisible(file.remove(f))

multiplyFiles= function(file1, file2, outfile, colsep="\t", overwrite=FALSE) {
  # Check args
  if ((!is.character(file1)) || (length(file1) != 1))
    stop("Argument 'file1' must be a file name.")
  if (!file.exists(file1))
    stop(paste("File '",file1,"' specified as 'file1' not found.",sep=""))
  if ((!is.character(file2)) || (length(file2) != 1))
    stop("Argument 'file2' must be a file name.")
  if (!file.exists(file2))
    stop(paste("File '",file2,"' specified as 'file2' not found.",sep=""))
  # Read tables using all character fields
  t1= read.table(file=file1, header=TRUE, sep=colsep, stringsAsFactors=FALSE)
  t2= read.table(file=file2, header=TRUE, sep=colsep, stringsAsFactors=FALSE)
  # Check compatibility of column headers
  if (!identical(names(t1),names(t2)))
    stop(paste("Column headers in '",file1,"' and '",file2,"' not identical.",sep=""))
  # Initialize result
  t3= t1
  # Try to multiply all columns
  for (i in 1:ncol(t1)) {
    if (is.numeric(t1[,i]) && is.numeric(t2[,i])) {
      t3[,i]= t1[,i] * t2[,i]
    } else if (is.character(t1[,i]) && is.character(t2[,i])) {
      if (!identical(t1[,i],t2[,i]))
        stop(paste("Column ",i," in file '",file1,"' and '",file2,"' is of type",
          " character but content is not identical.",sep=""))
    } else {
      stop(paste("Column ",i," in file '",file1,"' and '",file2,"' not compatible.",
        " Expecting either numbers or strings in all rows.",sep=""))
    }
  }
  # Write result
  if (file.exists(outfile) && (!overwrite))
    stop(paste("File '",outfile,"' already exists.",sep=""))
  write.table(x=t3, file=outfile, col.names=TRUE, row.names=FALSE, sep=colsep, quote=FALSE)
}


