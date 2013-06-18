#
#   ICF directory
#

L_ARCH = $(ARCH)

LIB_NAME = d-$(L_ARCH).a

BLAS = ../src/blas/$(LIB_NAME)

MINPACK2 = ../src/icf/$(LIB_NAME)

LIBS = $(MINPACK2) $(BLAS)

install: blas icf

blas:
	cd ../src/blas; make

icf:
	cd ../src/icf; make

CMEX = $(MATLABPATH)/bin/mex
CMEXFLAGS = -DARCH=$(ARCH) -largeArrayDims -g #-DMXDEBUG

# Path for mex.h

MINCLUDE = $(MATLABPATH)/extern/include
# MINCLUDE = /software/aix-4/com/packages/matlab-5.2/extern/include

matlab: lldl cgmex

lldl : lldl.c
	$(CMEX) $(CMEXFLAGS) -I$(MINCLUDE) lldl.c $(LIBS)

cgmex : cgmex.c
	$(CMEX) $(CMEXFLAGS) -I$(MINCLUDE) cgmex.c $(LIBS)


