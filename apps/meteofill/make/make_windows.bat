@echo off
setlocal

rem Name of the executable to be built
set BIN=../bin/meteofill

rem The directory of the c++ code collection containg a static library "libcpplib.a"
set cpplib=%ECHSE_GENERIC%/cpplib

rem Directory with C++ source code of the R package geostat
set geostatSrc=%ECHSE_TOOLS%/R/packages/geostat/src

rem Lists of the project-specific source and header files
set CPP=../src/meteofill.cpp ../src/meteofill_sub.cpp %geostatSrc%/geostat.cpp
set HDR=../src/meteofill_sub.h %geostatSrc%/geostat.h

rem Compile
g++ -ansi -iquote%geostatSrc% -iquote%cpplib% -L%cpplib% -lm -lstdc++ -pedantic -Wall -Wextra -ftrapping-math -ffast-math -fbounds-check -O3 -o %BIN% %CPP% -lcpplib

