\name{multiplyFiles}
\alias{multiplyFiles}
\title{Multiplication of text files of identical layout}
\usage{
multiplyFiles(file1, file2, outfile, colsep = "\\t", overwrite = FALSE)
}
\arguments{
  \item{file1}{First input file containing numeric data in
  one or more columns.  Column headers are mandatory.}

  \item{file2}{Second input file with identical layout as
  \code{file1}. In particular, column headers, the number
  of records, and the contents of any non-numeric columns
  must be identical.}

  \item{outfile}{Name/path of the result file.}

  \item{colsep}{Character to be used as a column separator
  in input and output files.}

  \item{overwrite}{Is it OK to replace an already existing
  output file?}
}
\value{
\code{NULL}.
}
\description{
The function multiplies corresponding values of the
\emph{numeric} columns of two text files. The result,
including non-numeric columns, is written to a result file.
}
\examples{
f= tempfile(c("in1","in2","out"))
t1= data.frame(city=c("Rome","Berlin","Oslo"),temperature=c(20,10,1))
t2= data.frame(city=c("Rome","Berlin","Oslo"),temperature=c(0.5,1,10))
write.table(t1, file=f[1], sep="\\t", row.names=FALSE, quote=FALSE)
write.table(t2, file=f[2], sep="\\t", row.names=FALSE, quote=FALSE)
multiplyFiles(file1=f[1], file2=f[2], outfile=f[3])
t3= read.table(file=f[3], header=TRUE)
print(t3)
invisible(file.remove(f))
}
\author{
David Kneis \email{david.kneis@uni-potsdam.de}
}

