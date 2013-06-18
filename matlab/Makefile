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

CMEX = $(MATLABDIR)/bin/mex
CMEXFLAGS = -DARCH=$(ARCH) -largeArrayDims -g #-DMXDEBUG

# Path for mex.h

MINCLUDE = $(MATLABDIR)/extern/include

matlab: lldl cgmex

lldl : lldl.c
	$(CMEX) $(CMEXFLAGS) -I$(MINCLUDE) lldl.c $(LIBS)

cgmex : cgmex.c
	$(CMEX) $(CMEXFLAGS) -I$(MINCLUDE) cgmex.c $(LIBS)


