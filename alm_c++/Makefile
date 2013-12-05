.SUFFIXES: .h .cpp

.cpp.o:
	${CXX} ${CXXFLAGS} -c $<

CXX = icpc -I /home/stsune/tadano/include
CXXFLAGS = -O2 -xSSE4.1 -mkl -vec-report1 -openmp -openmp-report1  -D_USE_EIGEN -D_VSL
LDFLAGS = -mkl

PROG = alm_c++

CXXSRC= alamode.cpp constraint.cpp error.cpp fcs.cpp files.cpp fitting.cpp input.cpp interaction.cpp \
	main.cpp memory.cpp symmetry.cpp system.cpp ewald.cpp timer.cpp writes.cpp patterndisp.cpp

OBJS= ${CXXSRC:.cpp=.o}

default: alm_c++

all: ${PROG}

alm_c++: ${OBJS}
	${CXX} ${LDFLAGS} -o $@ ${OBJS} ${LIBS}

clean:
	rm -f ${OBJS}

.PHONY: clean

depend:
	gcc -MM *.cpp > .depend

-include .depend
