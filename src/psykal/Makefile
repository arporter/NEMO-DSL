# Makefile for the original form of the tra_adv benchmark.

# Compiler and associated flags must be set through the following
# environment variables:
#
#   export F90=gfortran
#   export F90FLAGS=-O3
#   export OMPFLAGS=-fopenmp
#   export MPIF90=mpif90

all: tra_adv.exe

# Location of the dl_timer code
TIMER_DIR = ../shared/dl_timer

FORT_FLAGS = ${F90FLAGS} ${OMPFLAGS} -I${TIMER_DIR}/src

KERN_OBJ = zind_kern.o zwxy_kern.o

tra_adv.exe: timer_lib ${KERN_OBJ} psy.o tra_adv_alg.o
	${F90} ${OMPFLAGS} -o $@ tra_adv_alg.o psy.o ${KERN_OBJ} ${TIMER_DIR}/dl_timer_lib.a

%.o: %.f90
	${F90} ${FORT_FLAGS} -c $<
%.o: %.F90
	${F90} ${FORT_FLAGS} -c $<

timer_lib:
	make -C ${TIMER_DIR} sm_lib

clean:
	rm -f *.o
	rm -f *~
	rm -f *.mod
	rm -f tra_adv.exe
	rm -f output.dat
