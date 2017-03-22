# Makefile for the serial, psykal form of the tra_adv benchmark.

# Compiler and associated flags must be set through the following
# environment variables:
#
#   export F90=gfortran
#   export F90FLAGS="-O3 -openmp"

all: tra_adv.exe

# Location of the dl_timer code
TIMER_DIR = ../../shared/dl_timer

FORT_FLAGS = ${F90FLAGS} -I${TIMER_DIR}/src

tra_adv.exe: timer_lib tra_adv_kern.o tra_adv.o
	${F90} ${F90FLAGS} -o $@ tra_adv_kern.o tra_adv.o \
        ${TIMER_DIR}/dl_timer_lib.a

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
