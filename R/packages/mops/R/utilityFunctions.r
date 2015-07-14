#' Check existence of a directory
#'
#' The function checks whether a \emph{single} directory exists by evaluating
#' the result of a call to \code{\link{file.info}}.
#'
#' @param dirname A character string specifying the single path to be checked.
#'
#' @return A logical value.
#'
#' @seealso See \code{\link{file.exists}} for a similar test on files.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export

dir.exists = function(dirname) {
  stopifnot(length(dirname) == 1)
  x= file.info(dirname)$isdir
  if (is.na(x)) return(FALSE)
  return(x)
}


#' Move a single file to another directory
#'
#' The function moves a \emph{single} file to another directory using the
#' methods \code{\link{file.copy}} and \code{\link{file.remove}}.
#'
#' @param file The path of the file to be moved (character string).
#' @param dir The path of the target directory (character string).
#' @param overwrite A logical value. If \code{TRUE}, existing files in the
#'   target directory will silently be overwritten. Defaults to \code{FALSE}.
#'
#' @return \code{TRUE} on success and \code{FALSE} otherwise.
#'
#' @seealso See \code{\link{file.copy}} for copying.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export

file.move= function(file, dir, overwrite=FALSE) {
  stopifnot(length(file) == 1, length(dir) == 1, length(overwrite) == 1)
  if (!file.exists(file)) return(FALSE)
  if (!dir.exists(dir)) return(FALSE)
  if (!file.copy(from=file, to=dir, overwrite=overwrite)) return(FALSE)
  if (!file.remove(file)) return(FALSE)
  return(TRUE)
}

