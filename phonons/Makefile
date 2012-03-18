.SUFFIXES: .h .cpp

.cpp.o:
	${CXX} ${CXXFLAGS} -c $<

CXX = icpc
CXXFLAGS = -O2 -xSSE4.1 -vec-report1 
LDFLAGS = -mkl

PROG = phonons

CXXSRC= phonons.cpp error.cpp fcs_phonon.cpp parsephon.cpp dynamical.cpp  \
	main.cpp memory.cpp system.cpp timer.cpp write_phonons.cpp kpoint.cpp phonon_dos.cpp

OBJS= ${CXXSRC:.cpp=.o}

default: phonons

all: ${PROG}

phonons: ${OBJS}
	${CXX} ${LDFLAGS} -o $@ ${OBJS} ${LIBS}

clean:
	rm -f ${OBJS}

.PHONY: clean

depend:
	gcc -MM *.cpp > .depend

-include .depend
