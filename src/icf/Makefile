#
# ICF directory.
#

L_ARCH   = $(ARCH)
LIB_NAME = d-$(L_ARCH).a

OPTFLAGS = -O 

AR = ar -rcv

FILES = dicf.o   dpcg.o   dssyax.o   ihsort.o  srtdat2.o \
        dicfs.o  dsel.o   dsel2.o  dstrsol.o  insort.o

CFLAGS = $(OPTFLAGS) 
FFLAGS = $(OPTFLAGS)

icf: $(FILES) 
	$(AR) $(LIB_NAME) $(FILES)  

clean:
	- rm -f *.o

.c.o:
	$(CC) $(CFLAGS) -c $*.c
.f.o:
	$(FC) $(FFLAGS) -c $*.f


