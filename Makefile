#
# ICFS directory
#

L_ARCH   = $(ARCH)
LIB_NAME = d-$(L_ARCH).a

OPTFLAGS = -O

CFLAGS   = $(OPTFLAGS) 
FFLAGS   = $(OPTFLAGS)

# Libraries.

ICF    = src/icf/$(LIB_NAME)
BLAS   = src/blas/$(LIB_NAME)
TPROBS = src/tprobs/$(LIB_NAME)
UTILS   = src/utils/$(LIB_NAME)

LIBS = $(ICF) $(BLAS) $(UTILS)

install: libs exec

libs: icflib blas utils

icflib: 
	cd src/icf; make

blas: 
	cd src/blas; make

utils:
	cd src/utils; make

# Files for the MINPACK-2 incomplete Cholesky factorization.

exec : driver.o $(LIBS) 
	$(FC) $(FFLAGS) -o icf driver.o $(LIBS)

clean:
	cd src/icf;      make clean
	cd src/blas;     make clean
	cd src/utils;    make clean
	- rm -f *.o

.c.o:
	$(CC) $(CFLAGS) -c $*.c
.f.o:
	$(FC) $(FFLAGS) -c $*.f
