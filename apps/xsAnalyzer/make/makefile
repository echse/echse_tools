# Name of the executable to be built
BIN= ../bin/xsAnalyzer

# The directory of the c++ code collection containg a static library "libcpplib.a"
cpplib=$(ECHSE_GENERIC)/cpplib

# Lists of the project-specific source and header files
CPP= ../src/xsAnalyzer.cpp
HDR= 

# Compiler commands
CC     = /usr/bin/g++
CFLAGS = -ansi -iquote$(cpplib) -L$(cpplib) -lm -lstdc++ -pedantic -Wall -Wextra \
         -ftrapping-math -ffast-math -fbounds-check -O3

# The one and only target (default target)
all: $(OBJ) $(HDR) $(CPP)
	$(CC) $(CFLAGS) -o $(BIN) $(CPP) -lcpplib

