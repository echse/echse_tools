#' Code generator for the ECHSE
#'
#' Type \code{help(package="codegen")} to inspect the package description.
#'
#' @name codegen-package
#' @aliases codegen
#' @docType package
{}

#-------------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------------

# Expected column names in the input file
COLNAMES_EXPECTED= c("type","name")

# Long and short names of the various types of data involved in a object
# WARNING: The order of records in table 'TYPES' must not be changed as it
#          needs to conform with the order of the arguments in the constructor of
#          the template object group!
#          See file '?/echse_core/src/current/echse_coreClass_templateObjectGroup.h'
TYPES= data.frame(
  long=  c("paramFun","paramNum","inputExt","inputSim","stateScal","stateVect",
           "output","sharedParamFun","sharedParamNum"),
  desc=  c("Parameter function", "Scalar parameter", "External input", "Simulated input",
           "Scalar state variable", "Vector state variable", "Output variable",
           "Shared param. function", "Shared scalar param."),
  access_in_derivsScal= c(T,T,T,T,FALSE,T,T,T,T)
)

# Constants related to header files
FILE_CLASS_PREFIX= "AUTOechse_userClass_"
FILE_USER_PREFIX= "userCode_"
EXT_USER= "cpp"
EXT_HEADER= "h"
FILE_MAIN_INCL= "AUTOechse_includeFiles.h"
FILE_MAIN_FUNC= "AUTOechse_coreFunct_instantiateObjectGroups.h"

#-------------------------------------------------------------------------------
# Internal functions
#-------------------------------------------------------------------------------

# Test whether a string conforms with the naming rules for C++ identifiers
# (Names must start with a letter followed by letters, numbers, or "_".)
isValidName = function (s) {
  charset_1= letters
  charset_2= c(letters,"_",seq(from=0,to=9,by=1))
  if (nchar(s) < 1) return(F)
  s= tolower(s)
  if (!(substr(s,1,1) %in% charset_1)) return(F)
  if (nchar(s) > 1) {
    for (i in 2:nchar(s)) {
      if (!(substr(s,i,i) %in% charset_2)) return(F)
    }
  }
  return(T)
}

# Test whether a string is a reserved C++ keyword
isReservedWord = function (s) {
  reserved= c("and","and_eq","alignof","asm","auto","bitand","bitor","bool",
    "break","case","catch","char","char16_t","char32_t","class","compl","const",
    "constexpr","const_cast","continue","decltype","default","delete","double",
    "dynamic_cast","else","enum","explicittodo","export","externtodo","false",
    "float","for","friend","goto","if","inline","int","long","mutable",
    "namespace","new","noexcept","not","not_eq","null","nullptr","operator",
    "or","or_eq","private","protected","public","register","reinterpret_cast",
    "return","short","signed","sizeof","static","static_assert","static_cast",
    "struct","switch","template","this","thread_local","throw","true","try",
    "typedef","typeid","typename","union","using","unsigned","void","wchar_t",
    "virtual","volatile","while","xor","xor_eq")
  if (s %in% reserved) return(T)
  return(F)
}

# Append a record to a file
app= function(ofile,string) { write(x=string,file=ofile,ncolumns=1,append=T) }

# Initialize an output file
ofile_init = function (ofile,replace) {
  if ((!replace) && file.exists(ofile))
    stop(paste("Output file '",ofile,"' already exists.",sep=""))
  write(x=paste("// Created: ",date(),sep=""),file=ofile,ncolumns=1,append=F)
  app(ofile,"")
  app(ofile,"////////////////////////////////////////////////////////////////////////////////")
  app(ofile,"//////////////////// THIS IS A GENERATED FILE - DO NOT EDIT ////////////////////")
  app(ofile,"////////////////////////////////////////////////////////////////////////////////")
  app(ofile,"")
}

#-------------------------------------------------------------------------------
# Exported functions
#-------------------------------------------------------------------------------

#' Generates basic class code for an ECHSE-based model
#'
#' This function generates a set of C++ source files for use with the
#' Eco-Hydrological Simulation Environment (ECHSE). The source code is
#' generated based on a user-supplied specification of the model classes. For
#' each class, the specification includes the names and types of the class'
#' data members.
#'
#' @param files A \emph{named} vector of file names/paths. Each file is expected
#'   to contain the specification of a single model class. The class' name is
#'   determined by the name of the vector element. Class names must start with
#'   a letter, optionally followed by letters, digits, or underscores.
#'   Format and contents of these file(s) are described in detail below.
#' @param outdir Name/path of an \emph{existing} directory where all output
#'   files should be written to.
#' @param overwrite A logical value, specifying whether it is OK to overwrite
#'   existing files in the output directory given as \code{outdir}. Defaults to
#'   \code{FALSE}.
#'
#' @return NULL
#'
#' @note The files with the specification of model classes must be plain,
#'   TAB-separated text files with 2 mandatory columns 'type', and 'name'. Each
#'   file defines a separate class and each record declares a data member for
#'   the corresponding class. The name and type of the data member is taken from
#'   the 'name' and 'type' column, respectively. Member names must start with
#'   a letter, optionally followed by letters, digits, or underscores.
#'   Currently supported member types are 'stateScal' (for scalar
#'   state variables), 'stateVect' (for vector state variables), 'inputExt'
#'   (for external input variables, i.e. dynamic external forcing), 'inputSim'
#'   (for simulated input variable, i.e. dynamic forcing being the computed by
#'   another simulated object), 'paramNum' (for object-specific scalar numeric
#'   parameters), 'paramFun' (for object-specific parameter functions),
#'   'sharedParamNum' (for group-specific scalar numeric parameters),
#'   'sharedParamFun' (for group-specific parameter functions), and
#'   'output' (for output variables).
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#'
#' @examples
#' # Generation of a sample input files declaring two classes. The
#' # 1st class represents a linear reservoir, the 2nd class represents a link.
#' resv= data.frame(
#'   type= c("stateScal", "inputExt", "paramNum", "output"),
#'   name= c("volume",    "inflow",   "k",        "outflow")
#' )
#' f_resv= tempfile(pattern="classDecl_resv_",fileext=".txt")
#' write.table(x=resv, file=f_resv, sep="\t", row.names=FALSE, quote=FALSE)
#' link= data.frame(
#'   type= c("inputSim", "output"),
#'   name= c("inflow",   "outflow")
#' )
#' f_link= tempfile(pattern="classDecl_link_",fileext=".txt")
#' write.table(x=link, file=f_link, sep="\t", row.names=FALSE, quote=FALSE)
#' d= tempdir()
#' # Generation of C++ source files for an ECHSE-based model
#' generate(files=c(resv=f_resv, link=f_link), outdir=d)
#' print(paste("Outputs are in folder '",d,"'. Please remove the output",
#'   " manually after running this example.",sep=""))

generate= function(files, outdir, overwrite=FALSE) {
  # Check Args
  if (!is.character(files))
    stop("Argument 'files' must be a named vector of character strings.")
  if (is.null(names(files)) || (any(names(files) == "")))
    stop("Argument 'files' must be a named vector of character strings. Element name(s) missing.")
  if (any(duplicated(files)))
    stop("Elements in vector 'files' must be unique.")
  if (any(duplicated(names(files))))
    stop("Element names in vector 'files' must be unique.")
  for (f in files) {
    if (!file.exists(f))
      stop(paste("Input file '",f,"' does not exist.",sep=""))
  }
  if ((!is.finite(overwrite)) || (!is.logical(overwrite)) || (length(overwrite)!=1))
    stop("Argument 'overwrite' must be TRUE or FALSE.")
  if ((!is.character(outdir)) || (length(outdir) != 1))
    stop("Argument 'outdir' must be a single path name.")
  if (!file.exists(outdir))
    stop(paste("Output directory '",outdir,"' does not exist.",sep=""))
  if (is.na(file.info(outdir)$isdir))
    stop(paste("Path '",outdir,"' does not point to a directory.",sep=""))

  #-----------------------------------------------------------------------------
  # Read input table(s) with class declaration and put into a single table
  for (i in 1:length(files)) {
    print(paste("Reading declaration of class '",names(files[i]),"' from '",
      basename(files[i]),"...",sep=""))
    t= read.table(file=files[i],header=T,sep="\t",quote="",colClasses=c("character"))
    names(t)= tolower(names(t))
    if (!all(COLNAMES_EXPECTED %in% names(t))) {
      stop(paste("Missing column(s) in input file '",files[i],"'. Expected: '",
        paste(COLNAMES_EXPECTED,collapse="', "),"'.",sep=""))
    }
    t$class= names(files[i])  # Set class name from name of element of vector 'files'
    if (i == 1) {
      tab= t
    } else {
      if (!identical(names(t),names(tab)))
        stop(paste("Column names in '",files[i],"' not consistent with names in '",
          files[1],"'.",sep=""))
      tab= rbind(tab, t)
    }
  }

  # Check that names of classes and features are valid C++ identifier names
  print("Checking class and member names...")
  for (i in 1:nrow(tab)) {
    if (!isValidName(tab$class[i])) stop(paste("String '",tab$class[i],
      "' cannot be used as a class name as it is not a valid C++ identifier.",
      " Identifier names must start with a letter followed by letter(s),",
      " number(s) or unserscore(s).",sep=""))
    if (isReservedWord(tab$class[i])) stop(paste("String '",tab$class[i],
      "' cannot be used as a class name as it is a reserved keyword in C++.",sep=""))
    if (!isValidName(tab$name[i])) stop(paste("String '",tab$name[i],
      "' cannot be used as a feature name in class '",tab$class[i],
      "' as it is not a valid C++ identifier.",
      " Identifier names must start with a letter followed by letter(s),",
      " number(s) or underscore(s).",sep=""))
    if (isReservedWord(tab$name[i])) stop(paste("String '",tab$name[i],
      "' cannot be used as a feature name in class '",tab$class[i],
      "' as it is reserved keyword in C++.",sep=""))
  }
  # Check that only valid feature types are present
  print("Checking member types...")
  for (i in 1:nrow(tab)) {
    if (!(tab$type[i] %in% TYPES$long)) stop(paste("Member type '",tab$type[i],
      "' used in class '",tab$class[i],"' is unknown. Type must be one of: '",
      paste(TYPES$long,collapse="', "),"'.",sep=""))
  }
  # Check that feature names within a class are unique
  print("Checking member names for uniqueness...")
  for (cls in unique(tab$class)) {
    names= unique(tab[tab$class==cls, "name"])
    if (length(unique(names)) != nrow(tab[tab$class==cls,])) {
      stop(paste("Duplicate feature names in class '",cls,
        "' detected. Names must be unique within the class.",sep=""))
    }
  }

  print("Creating output files...")

  #-----------------------------------------------------------------------------
  # Creation of header that collects (includes) the headers of all object classes
  ofile= paste(outdir,"/",FILE_MAIN_INCL,sep="")
  ofile_init(ofile,overwrite)
  app(ofile,"// This file is to bundle various include files into a single one")
  app(ofile,"")
  app(ofile,"// The file below provides the function to instantiate the objectGroups")
  app(ofile,paste("#include \"",FILE_MAIN_FUNC,"\"",sep=""))
  app(ofile,"")
  app(ofile,"// The files below provide declarations of user-defined object classes")
  for (cls in unique(tab$class)) {
    app(ofile,paste("#include \"",FILE_CLASS_PREFIX,cls,".",EXT_HEADER,"\"",sep=""))
  }
  app(ofile,"")

  #-----------------------------------------------------------------------------
  # Creation of header files for object classes
  for (cls in unique(tab$class)) {
    ofile= paste(outdir,"/",FILE_CLASS_PREFIX,cls,".",EXT_HEADER,sep="")
    ofile_init(ofile,overwrite)
    app(ofile,paste("// This file provides the declaration of class '",cls,"'.",sep=""))
    app(ofile,"//")
    app(ofile,"// To make this class usable, this header files must be supplemented")
    app(ofile,"// by the body code of the class methods. The code is assumed to")
    app(ofile,"// reside in the following files:")
    app(ofile,paste("//   <some_directory>/",FILE_USER_PREFIX,cls,"_aux.",EXT_USER,sep=""))
    app(ofile,paste("//   <some_directory>/",FILE_USER_PREFIX,cls,"_simulate.",EXT_USER,sep=""))
    app(ofile,paste("//   <some_directory>/",FILE_USER_PREFIX,cls,"_derivsScal.",EXT_USER,sep=""))
    app(ofile,"// Note that the last file may be empty, if the 'derivsScal' method")
    app(ofile,"// of this class is not used. The first file may contain definitions")
    app(ofile,"// (of functions, for example) that should be available in the two methods.")
    app(ofile,"// This file may also be empty.")
    app(ofile,"")
    app(ofile,paste("#ifndef ",toupper(FILE_CLASS_PREFIX),toupper(cls),"_",toupper(EXT_HEADER),sep=""))
    app(ofile,paste("#define ",toupper(FILE_CLASS_PREFIX),toupper(cls),"_",toupper(EXT_HEADER),sep=""))
    app(ofile,"")
    app(ofile,"// Include standard headers")
    app(ofile,"#include <cmath>    // Math functions")
    app(ofile,"#include <sstream>  // For generating error messages")
    app(ofile,"#include \"echse_coreClass_abstractObject.h\" // Abstract base class")
    app(ofile,"#include \"echse_coreFunct_solveODE.h\"       // Access to ODE solver")
    app(ofile,"")
    app(ofile,"using namespace std;")
    app(ofile,"")
    app(ofile,paste("class object_",cls,": public abstractObject {",sep=""))
    app(ofile,"  public:")
    app(ofile,paste("    object_",cls,"() { };   // Constructor - Definition not required",sep=""))
    app(ofile,paste("    ~object_",cls,"() { };  // Destructor - Definition not required",sep=""))
    app(ofile,"    // Virtual methods inherited from abstract base class")
    app(ofile,"    void simulate(const unsigned int delta_t);")
    app(ofile,"    void derivsScal(const double t, const vector<double> &u,")
    app(ofile,"      vector<double> &dudt, const unsigned int delta_t);")
    app(ofile,"};")
    app(ofile,"")
    app(ofile,"// ----- Frames of method implementations follow -----")
    app(ofile,"")
    app(ofile,"// Import auxiliary definitions to be used in 'simulate' and 'derivsScal'")
    app(ofile,paste("#include \"",FILE_USER_PREFIX,cls,"_aux.",EXT_USER,"\"",sep=""))
    app(ofile,"")
    app(ofile,"// The 'simulate' method")
    app(ofile,paste("namespace object_",cls,"_simulate {",sep=""))
    for (it in 1:nrow(TYPES)) {
      # app(ofile,paste("  // ",TYPES$desc[it]," ('",TYPES$long[it],"')",sep=""))
      n= tab$name[tab$class==cls & tab$type==TYPES$long[it]]
      if (length(n) > 0) {
        for (i in 1:length(n)) {
          app(ofile,paste("  static const T_index_",TYPES$long[it]," ",n[i],"= {",i-1,"};", sep=""))
        }
      }
    }
    app(ofile,"}")
    app(ofile,paste("void object_",cls,"::simulate(const unsigned int delta_t) {",sep=""))
    app(ofile,paste("  using namespace object_",cls,"_simulate;",sep=""))
    app(ofile,"  try {")
    app(ofile,paste("    #include \"",FILE_USER_PREFIX,cls,"_simulate.",EXT_USER,"\"",sep=""))
    app(ofile,"  } catch (except) {")
    app(ofile,"    stringstream errmsg;")
    app(ofile,"    errmsg << \"Simulation failed for time step of length \" << delta_t << \".\";")
    app(ofile,paste("    except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);",sep=""))
    app(ofile,"    throw(e);")
    app(ofile,"  }")
    app(ofile,"}")
    app(ofile,"")
    app(ofile,"// The 'derivsScal' method")
    app(ofile,paste("namespace object_",cls,"_derivsScal {",sep=""))
    for (it in 1:nrow(TYPES)) {
      # app(ofile,paste("  // ",TYPES$desc[it]," ('",TYPES$long[it],"')",sep=""))
      n= tab$name[tab$class==cls & tab$type==TYPES$long[it]]
      if (length(n) > 0) {
        for (i in 1:length(n)) {
          if (TYPES$access_in_derivsScal[it]) {
            app(ofile,paste("  static const T_index_",TYPES$long[it]," ",n[i],"= {",i-1,"};", sep=""))
          } else {
            app(ofile,paste("  static const unsigned int INDEX_",n[i],"= ",i-1,"; // for use with vectors u & dudt", sep=""))
          }
        }
      }
    }
    app(ofile,"}")
    app(ofile,paste("void object_",cls,"::derivsScal(const double t,",sep=""))
    app(ofile,"  const vector<double> &u, vector<double> &dudt, const unsigned int delta_t) {")
    app(ofile,paste("  using namespace object_",cls,"_derivsScal;",sep=""))
    app(ofile,"  try {")
    app(ofile,paste("    #include \"",FILE_USER_PREFIX,cls,"_derivsScal.",EXT_USER,"\"",sep=""))
    app(ofile,"  } catch (except) {")
    app(ofile,"    stringstream errmsg, data;")
    app(ofile,"    for (unsigned int i=0; i<u.size(); i++) data << u[i] << \" \";")
    app(ofile,paste("    errmsg << \"Calculation of derivatives failed (Time: \"",
      " << t << \" Number\" <<\n      \" of deriv.s: \" << dudt.size() << \"",
      " Value(s) of state(s): \" << data.str() <<\n      \" Length of time step: \"",
      " << delta_t << \").\";",sep=""))
    app(ofile,paste("    except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);",sep=""))
    app(ofile,"    throw(e);")
    app(ofile,"  }")
    app(ofile,"}")
    app(ofile,"")
    app(ofile,"#endif")
    app(ofile,"")
  }

  #-----------------------------------------------------------------------------
  # Creation of object group classes from template class
  ofile= paste(outdir,"/",FILE_MAIN_FUNC,sep="")
  ofile_init(ofile,overwrite)
  app(ofile,"// The function defined below returns an array of the objectGroups.")
  app(ofile,"//")
  app(ofile,"// Details: For each objectGroup, a single instance is created and")
  app(ofile,"// pointers to these instances are stored in the output array.")
  app(ofile,"// Note that the objectGroup instances must be declared as STATIC")
  app(ofile,"// for the pointers to be valid outside the function.")
  app(ofile,"")
  app(ofile,"// Base headers")
  app(ofile,"#include <string>")
  app(ofile,"#include <vector>")
  app(ofile,"")
  app(ofile,"// Constant headers")
  app(ofile,"#include \"echse_coreClass_abstractObjectGroup.h\"")
  app(ofile,"#include \"echse_coreClass_templateObjectGroup.h\"")
  app(ofile,"")
  app(ofile,"// Generated headers")
  for (cls in unique(tab$class)) {
    app(ofile,paste("#include \"",FILE_CLASS_PREFIX,cls,".",EXT_HEADER,"\"",sep=""))
  }
  app(ofile,"")
  app(ofile,"void instantiateObjectGroups(")
  app(ofile,"  vector<abstractObjectGroup*> &objectGroups")
  app(ofile,") {")
  app(ofile,"")
  for (it in TYPES$long) {
    app(ofile,paste("  vector<string> ",it,";",sep=""))
  }
  app(ofile,"")
  app(ofile,"  // PART 1: Instantiation of objectGroups from template class")
  app(ofile,"")
  for (cls in unique(tab$class)) {
    app(ofile,paste("  // OBJECT GROUP '",cls,"':",sep=""))
    for (it in TYPES$long) {
      app(ofile,paste("  ",it," = {",sep=""))
      n= tab$name[tab$class==cls & tab$type==it]
      if (length(n) > 0) {
        for (i in 1:length(n)) {
          if (i < length(n)) {
            app(ofile,paste("     \"",n[i],"\",",sep=""))
          } else {
            app(ofile,paste("     \"",n[i],"\"",sep=""))
          }
        }
      }
      app(ofile,"  };")
    }
    app(ofile,paste("  static templateObjectGroup<object_",cls,
      "> objectGroup_",cls,"(\"",cls,"\",",sep=""))
    for (it in TYPES$long) {
      if (it != TYPES$long[nrow(TYPES)]) {
        app(ofile,paste("    ",it,",",sep=""))
      } else {
        app(ofile,paste("    ",it,sep=""))
      }
    }
    app(ofile,"  );")
    app(ofile," ")
  }
  app(ofile,"")
  app(ofile,"  // PART 2: Set array of pointers to static objectGroups")
  app(ofile,"")
  mfs= unique(tab$class)
  app(ofile,paste("  objectGroups.resize(",length(mfs),");",sep=""))
  for (i in 1:length(mfs)) {
    app(ofile,paste("  objectGroups[",(i-1),"]= &objectGroup_",mfs[i],";",sep=""))
  }
  app(ofile,"} // End of function\n")

  #-----------------------------------------------------------------------------
  return(invisible(NULL))
}

