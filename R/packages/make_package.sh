#!/bin/bash

# This is a generic script to build, check, and install an R package.
#
# It is assumed that the package's documentation is contained in the source code
# for processing with roxygen.

## BEGIN SETTINGS ##############################################################

# A folder where R is searching for packages --> Adjust name of folder and R version
installDir=$HOME"/R/i686-pc-linux-gnu-library/3.1"

## END SETTINGS ################################################################


# Get name and version of package to be build
if [ -z $1 -o -z $2 ]; then
  echo Usage: $0 name_of_package version_number
  exit 1
fi
package=$1
version=$2

# Set name of folder containing the built package
build="$ECHSE_TOOLS/R/packages/$package.Rcheck/$package"

# Create documentation generation script to be processed by R
echo "Creating documentation files..."
f=$(mktemp)
cat > $f << EOL
library('roxygen2')
roxygenize('${package}')
EOL
cat $f
# Let R create the documentation from in-source comments
R --vanilla --silent --slave --file=$f
if [ $? -ne 0 ]; then
  rm $f
  echo Failed to create the documentation files.
  exit 1
else
  rm $f
fi

# Build and check the package
echo "Building package..."
R CMD build $package
if [ $? -ne 0 ]; then
  echo Failed to build the package.
  exit 1  
fi
echo "Checking package..."
R CMD check $package"_"$version".tar.gz"
if [ $? -ne 0 ]; then
  echo Package check returned non-zero status.
  exit 1  
fi

# Move package to R's install dir
echo "Installing package..."
if [ -d "$installDir/$package" ]; then
  rm -r "$installDir/$package"
fi
mv $build $installDir
if [ $? -ne 0 ]; then
  echo Failed to move the built package to $installDir.
  exit 1  
fi

# Normal exit
exit 0

