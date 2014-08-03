#
# ICFS directory
#
SHELL = /bin/bash

L_ARCH   = $(ARCH)
LIB_NAME = liblldl.dylib
LIBUTILS_NAME = liblldl_utils.dylib

OPTFLAGS = -O

CFLAGS   = $(OPTFLAGS) -fPIC
FFLAGS   = $(OPTFLAGS) -fPIC

# Libraries.

ICF    = src/icf/$(LIB_NAME)
TPROBS = src/tprobs/$(LIB_NAME)
UTILS   = src/utils/$(LIBUTILS_NAME)

LIBS = $(ICF) $(UTILS)

PREFIX ?= $(PWD)

install: libs exec
	[[ ! -d $(PREFIX)/bin ]] && mkdir -p $(PREFIX)/bin || true
	[[ ! -d $(PREFIX)/lib ]] && mkdir -p $(PREFIX)/lib || true
	cp icf $(PREFIX)/bin
	cp src/icf/liblldl.dylib src/utils/liblldl_utils.dylib $(PREFIX)/lib

libs: icflib utils

icflib: 
	cd src/icf; CFLAGS="$(CFLAGS)" FFLAGS="$(FFLAGS)" make

utils:
	cd src/utils; CFLAGS="$(CFLAGS)" FFLAGS="$(FFLAGS)" make

# Files for the MINPACK-2 incomplete Cholesky factorization.

exec : driver.o $(LIBS) 
	$(FC) $(FFLAGS) -o icf driver.o $(LIBS)

clean:
	cd src/icf;      make clean
	cd src/utils;    make clean
	- rm -f *.o

.c.o:
	$(CC) $(CFLAGS) -c $*.c
.f.o:
	$(FC) $(FFLAGS) -c $*.f
