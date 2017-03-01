# Makefile for the PSyKAl form of the tra_adv benchmark.

# Compiler and associated flags must be set through the following
# environment variables:
#
#   export F90=gfortran
#   export F90FLAGS=-O3

all: tra_adv.exe

# Location of the dl_timer code
TIMER_DIR = ../shared/dl_timer
INF_DIR = ../shared/dl_esm_inf/finite_difference

# Append the locations of the various library modules to the
# list of compiler flags
FORT_FLAGS = ${F90FLAGS} -I${TIMER_DIR}/src -I${INF_DIR}/src

KERN_OBJ = zind_kern.o zwxy_kern.o zslpxy_kern.o zslpxy_update_kern.o zwxy2_kern.o \
	   mydomain_update_kern.o zwx_kern.o zslpx_kern.o zslpx_update_kern.o \
	   zwx2_kern.o mydomain_kern.o

tra_adv.exe: timer_lib inf_lib ${KERN_OBJ} psy.o tra_adv_alg.o
	${F90} -o $@ tra_adv_alg.o psy.o ${KERN_OBJ} \
        ${TIMER_DIR}/dl_timer_lib.a ${INF_DIR}/src/dl_esm_inf_fd.a

%.o: %.f90
	${F90} ${FORT_FLAGS} -c $<
%.o: %.F90
	${F90} ${FORT_FLAGS} -c $<

timer_lib:
	make -C ${TIMER_DIR} sm_lib

inf_lib:
	make -C ${INF_DIR}

clean:
	rm -f *.o
	rm -f *~
	rm -f *.mod
	rm -f tra_adv.exe
	rm -f output.dat
