.SUFFIXES: .h .cpp

.cpp.o:
	${CXX} ${CXXFLAGS} -c $<

CXX = icpc
CXXFLAGS = -O3 -xSSE4.1 -vec-report1 -C -debug all
LDFLAGS = -mkl

PROG = alm_c++

CXXSRC= alamode.cpp error.cpp fcs.cpp files.cpp fitting.cpp input.cpp interaction.cpp \
	main.cpp memory.cpp symmetry.cpp system.cpp timer.cpp

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
