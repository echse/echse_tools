\name{run_model}
\alias{run_model}
\title{Run an external model}
\usage{
run_model(model_path, model_args, stopIfNonZero = FALSE)
}
\arguments{
  \item{model_path}{Path of the model executable. Note that
  the executable is executed in R's current working
  directory, i.e. the one returned by \code{\link{getwd}}.
  This is true even if the executable resides somewhere
  else.}

  \item{model_args}{A \emph{named} vector of command line
  arguments to be passed to the model. Note that, for named
  elements, name and value are concatenated by '=',
  automatically. Moreover, all values are enclosed in shell
  quotes (using \code{\link{shQuote}}). For example, a
  vector \code{c(a='1', '-x', b='2')} is expanded to the
  string 'a='1' '-x' b='2''.}

  \item{stopIfNonZero}{If \code{TRUE} the function calls
  \code{stop} if the model executable returns a non-zero
  exit code. If \code{FALSE}, only a warning is generated.}
}
\value{
The return code of the model executable.
}
\description{
The function runs a model (being an executable file) with
the supplied set of command line arguments and captures the
return code.
}
\author{
David Kneis \email{david.kneis@uni-potsdam.de}
}

