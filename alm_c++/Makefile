.SUFFIXES: .h .cpp

.cpp.o:
	${CXX} ${CXXFLAGS} -c $<

CXX = icpc -I /home/stsune/tadano/include
CXXFLAGS = -O2 -xSSE4.1 -vec-report1 -openmp -openmp-report1 
LDFLAGS = -mkl

PROG = alm_c++

CXXSRC= alamode.cpp constraint.cpp error.cpp fcs.cpp files.cpp fitting.cpp input.cpp interaction.cpp \
	main.cpp memory.cpp symmetry.cpp system.cpp timer.cpp writes.cpp

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
