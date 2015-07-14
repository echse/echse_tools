
This is a software tool which can be used to fill gaps in multi-location
time series of meteorological variables by spatial interpolation.

(C) 2013, David Kneis, david.kneis@uni-potsdam.de

NOTES:
------

This software is an application written in C++ and must be compiled before it can
be used. The recommended compiler is 'g++' from the GNU compiler collection.
The 'make' directory contains a makefile for LINUX systems and an equivalent
batch file for compilation on WINDOWS systems using g++.

IMPORTANT: The compilation requires access to a static library 'libcpplib.a'
which is distributed with the ECHSE modeling software. For successful compilation
it is essential that this library is (1) up-to-date and (2) compiled with the
same compiler/version as this program. If in doubt, please re-compile the
library before compiling this program.

See the file 'ChangeLog.txt' for information on the latest modifications.

