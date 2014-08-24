.SUFFIXES: .h .cpp
#-----------------------------------------------
# Makefile for 'tools' directory
# Please modify the variables properly.
# We recommend to use Intel c++ compiler.
#-----------------------------------------------

#Intel compiler
#CXX = icpc

# OSX
CXX = c++

CXXFLAGS = -O2 
INCLUDE = -I../include 

#----------------------------------------------
# General rules
#----------------------------------------------

.cpp.o:
	${CXX} ${CXXFLAGS} ${INCLUDE} -c $<

PROG = analyze_phonons

CXXSRC = analyze_phonons.cpp

OBJS= ${CXXSRC:.cpp=.o}

default: all

all: ${PROG}

analyze_phonons: ${OBJS}
	${CXX} ${CXXFLAGS} -o $@ ${OBJS} 

clean:
	rm -f ${OBJS}

.PHONY: clean

depend:
	gcc -MM *.cpp > .depend

-include .depend
