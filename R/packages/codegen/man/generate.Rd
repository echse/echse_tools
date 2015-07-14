\name{generate}
\alias{generate}
\title{Generates basic class code for an ECHSE-based model}
\usage{
generate(files, outdir, overwrite = FALSE)
}
\arguments{
  \item{files}{A \emph{named} vector of file names/paths.
  Each file is expected to contain the specification of a
  single model class. The class' name is determined by the
  name of the vector element. Class names must start with a
  letter, optionally followed by letters, digits, or
  underscores.  Format and contents of these file(s) are
  described in detail below.}

  \item{outdir}{Name/path of an \emph{existing} directory
  where all output files should be written to.}

  \item{overwrite}{A logical value, specifying whether it
  is OK to overwrite existing files in the output directory
  given as \code{outdir}. Defaults to \code{FALSE}.}
}
\description{
This function generates a set of C++ source files for use
with the Eco-Hydrological Simulation Environment (ECHSE).
The source code is generated based on a user-supplied
specification of the model classes. For each class, the
specification includes the names and types of the class'
data members.
}
\note{
The files with the specification of model classes must be
plain, TAB-separated text files with 2 mandatory columns
'type', and 'name'. Each file defines a separate class and
each record declares a data member for the corresponding
class. The name and type of the data member is taken from
the 'name' and 'type' column, respectively. Member names
must start with a letter, optionally followed by letters,
digits, or underscores.  Currently supported member types
are 'stateScal' (for scalar state variables), 'stateVect'
(for vector state variables), 'inputExt' (for external
input variables, i.e. dynamic external forcing), 'inputSim'
(for simulated input variable, i.e. dynamic forcing being
the computed by another simulated object), 'paramNum' (for
object-specific scalar numeric parameters), 'paramFun' (for
object-specific parameter functions), 'sharedParamNum' (for
group-specific scalar numeric parameters), 'sharedParamFun'
(for group-specific parameter functions), and 'output' (for
output variables).
}
\examples{
# Generation of a sample input files declaring two classes. The
# 1st class represents a linear reservoir, the 2nd class represents a link.
resv= data.frame(
  type= c("stateScal", "inputExt", "paramNum", "output"),
  name= c("volume",    "inflow",   "k",        "outflow")
)
f_resv= tempfile(pattern="classDecl_resv_",fileext=".txt")
write.table(x=resv, file=f_resv, sep="\\t", row.names=FALSE, quote=FALSE)
link= data.frame(
  type= c("inputSim", "output"),
  name= c("inflow",   "outflow")
)
f_link= tempfile(pattern="classDecl_link_",fileext=".txt")
write.table(x=link, file=f_link, sep="\\t", row.names=FALSE, quote=FALSE)
d= tempdir()
# Generation of C++ source files for an ECHSE-based model
generate(files=c(resv=f_resv, link=f_link), outdir=d)
print(paste("Outputs are in folder '",d,"'. Please remove the output",
  " manually after running this example.",sep=""))
}
\author{
David Kneis \email{david.kneis@uni-potsdam.de}
}

