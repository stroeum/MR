#This makefile create the executable file FractalModel
CC    = CC
CXX   = mpicxx
RUN   = mpiexec
MAKE  = MAKE

#FLAGS = -Os -g -Wall -std=gnu99
FLAGS = -O3 -ffast-math -Wall
CXXFLAGS = $(FLAGS)

all: all-redirect

PROGRAMS = test

SOURCES = MathFunctions.cpp AtmScaling.cpp MpiInput.cpp Matrix.cpp ResGrid.cpp SizeGrid.cpp SizeDomain.cpp Sources.cpp ScaledFields.cpp MatrixFunctions.cpp Cloud.cpp ConductivityProfile.cpp SorSolution.cpp Lax.cpp Input.cpp MpiFunctions.cpp IOFunctions.cpp BoundaryConditions.cpp

OBJECTS = MathFunctions.o AtmScaling.o MpiInput.o Matrix.o ResGrid.o SizeGrid.o SizeDomain.o Sources.o ScaledFields.o MatrixFunctions.o Cloud.o ConductivityProfile.o SorSolution.o Lax.o Input.o MpiFunctions.o IOFunctions.o BoundaryConditions.o

PGM_SOURCES = ${SOURCES} main.cpp

${PROGRAMS}: main.o ${OBJECTS}
	${CXX} $(CXXFLAGS)  -lm -o ${PROGRAMS} main.o ${OBJECTS}

all-redirect: ${PROGRAMS}
	#Automatically execute the compiled code
	${RUN} -n 8 ./${PROGRAMS}
	#rm -rf *.o ${PROGRAMS}	results/*.dat results/*.mat results/*.avi

clear:
	rm -rf results/phi2d1*.dat
	rm -rf results/phi2d2*.dat
	rm -rf results/phi2d3*.dat
	rm -rf results/phi2d4*.dat
	rm -rf results/phi2d5*.dat
	rm -rf results/phi2d*.dat
	rm -rf results/rhos2d1*.dat
	rm -rf results/rhos2d2*.dat
	rm -rf results/rhos2d3*.dat
	rm -rf results/rhos2d4*.dat
	rm -rf results/rhos2d5*.dat
	rm -rf results/rhos2d*.dat
	rm -rf results/rho2d1*.dat
	rm -rf results/rho2d2*.dat
	rm -rf results/rho2d3*.dat
	rm -rf results/rho2d4*.dat
	rm -rf results/rho2d5*.dat
	rm -rf results/rho2d*.dat
	rm -rf results/Er2d1*.dat
	rm -rf results/Er2d2*.dat
	rm -rf results/Er2d3*.dat
	rm -rf results/Er2d4*.dat
	rm -rf results/Er2d5*.dat
	rm -rf results/Er2d*.dat
	rm -rf results/Ez2d1*.dat
	rm -rf results/Ez2d2*.dat
	rm -rf results/Ez2d3*.dat
	rm -rf results/Ez2d4*.dat
	rm -rf results/Ez2d5*.dat
	rm -rf results/Ez2d*.dat
	rm -rf *.o ${PROGRAMS} *~ PI* results/*.dat results/*.mat results/*.avi

.SUFFIXES: .o .c .f .cc .f90
.c:
	${CC} ${FLAGS} -o $* $< 
.c.o:
	${CC} ${FLAGS} -c $<
.cc:
	${CXX} ${FLAGS} -o $* $<  

sources:
	@echo ${PGM_SOURCES}