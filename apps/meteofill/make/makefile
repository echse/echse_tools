# Name of the executable to be built
BIN= ../bin/meteofill

# The directory of the c++ code collection containg a static library "libcpplib.a"
cpplib=$(ECHSE_GENERIC)/cpplib

# Directory with C++ source code of the R package geostat
geostatSrc=$(ECHSE_TOOLS)/R/packages/geostat/src

# Lists of the project-specific source and header files
CPP= ../src/meteofill.cpp ../src/meteofill_sub.cpp $(geostatSrc)/geostat.cpp
HDR= ../src/meteofill_sub.h $(geostatSrc)/geostat.h

# Compiler commands
CC     = /usr/bin/g++
CFLAGS = -ansi -iquote$(geostatSrc) -iquote$(cpplib) -L$(cpplib) -lm -lstdc++ -pedantic -Wall -Wextra \
         -ftrapping-math -ffast-math -fbounds-check -O3

# The one and only target (default target)
all: $(OBJ) $(HDR) $(CPP)
	$(CC) $(CFLAGS) -o $(BIN) $(CPP) -lcpplib

