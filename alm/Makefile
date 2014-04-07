.SUFFIXES: .h .cpp
#-----------------------------------------------
# Makefile for the program 'ALM'.
# Please modify the variables properly.
# We recommend to use Intel c++ compiler.
#-----------------------------------------------

CXX = icpc 
CXXFLAGS = -O2 -xHOST -openmp 
#CXXFLAGS = -O2 -xHOST -openmp -D_USE_EIGEN
INCLUDE = -I../include -I/home/stsune/tadano/include

CXXL = ${CXX}
LDFLAGS = -mkl

LAPACK = 
LIBS = ${LAPACK}

#-----------------------------------------------
# General rules 
#-----------------------------------------------

.cpp.o:
	${CXX} ${CXXFLAGS} ${INCLUDE} -c $<

PROG = alm

CXXSRC= alamode.cpp constraint.cpp error.cpp fcs.cpp files.cpp \
	fitting.cpp input.cpp interaction.cpp main.cpp memory.cpp \
	patterndisp.cpp symmetry.cpp system.cpp timer.cpp writes.cpp 

OBJS= ${CXXSRC:.cpp=.o}

default: alm

all: ${PROG}

alm: ${OBJS}
	${CXXL} ${LDFLAGS} -o $@ ${OBJS} ${LIBS}

clean:
	rm -f ${OBJS}

.PHONY: clean

depend:
	gcc -MM *.cpp > .depend

-include .depend
