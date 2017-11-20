.SUFFIXES: .h .cpp
#-----------------------------------------------
# Makefile for 'tools' directory
# Please modify the variables properly.
# We recommend to use Intel c++ compiler.
#-----------------------------------------------

#Intel compiler
#CXX = icpc

# OSX
CXX = g++-7

CXXFLAGS = -O2 

# Add the directory containing boost header files here
INCLUDE = -I../include -I$(HOME)/include

#----------------------------------------------
# General rules
#----------------------------------------------

.cpp.o:
	${CXX} ${CXXFLAGS} ${INCLUDE} -c $<

PROG = analyze_phonons qe2alm dfc2

CXXSRC = analyze_phonons.cpp qe2alm.cpp dfc2.cpp

OBJS= ${CXXSRC:.cpp=.o}

default: all

all: ${PROG}

analyze_phonons: analyze_phonons.o
	${CXX} ${CXXFLAGS} -o $@ $<

qe2alm: qe2alm.o
	${CXX} ${CXXFLAGS} -o $@ $<

dfc2: dfc2.o
	${CXX} ${CXXFLAGS} -o $@ $<

clean:
	rm -f ${OBJS}

.PHONY: clean

depend:
	gcc -MM *.cpp > .depend

-include .depend
