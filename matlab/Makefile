#
#   ICF directory
#

L_ARCH = $(ARCH)

LIB_NAME = liblldl.dylib

MINPACK2 = ../src/icf/$(LIB_NAME)

LIBS = $(MINPACK2)

install: icf

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


