#!/bin/bash -i

datafile=$1

thisDir="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rscript=$thisDir"/tsplot.r"
errfile=$datafile".tsplot_log"

cd $(dirname "$datafile")
R --vanilla --file=$rscript --args $datafile 2>$errfile

if [ ! -s $errfile ]; then
  rm $errfile
fi

