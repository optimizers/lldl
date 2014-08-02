#
# ICF directory.
#

L_ARCH   = $(ARCH)
LIB_NAME = liblldl.dylib

AR = $(FC) -shared -Wl,-framework -Wl,Accelerate -Wl,-install_name -Wl,$(PREFIX)/lib/$(LIB_NAME) -o

FILES = dicf.o   dpcg.o   dssyax.o   ihsort.o  srtdat2.o \
        dicfs.o  dsel.o   dsel2.o  dstrsol.o  insort.o

icf: $(FILES)
	$(AR) $(LIB_NAME) $(FILES)

clean:
	- rm -f *.o

.c.o:
	$(CC) $(CFLAGS) -c $*.c
.f.o:
	$(FC) $(FFLAGS) -c $*.f


