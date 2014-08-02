#
# UTILS library.
#

L_ARCH   = $(ARCH)
LIB_NAME = liblldl_utils.dylib

CPPFLAGS  = -E -C -D$(ARCH) -fPIC

AR = $(FC) -shared -Wl,-framework -Wl,Accelerate -Wl,-install_name -Wl,$(PREFIX)/lib/$(LIB_NAME) -o

FILES = cputime_$(ARCH).o  dpmeps.o  mmio.o

src: $(FILES)
	$(AR) $(LIB_NAME) $(FILES)

clean:
	- rm -f *.o

.SUFFIXES: .c .f .F

.c.o:
	$(CC) $(CFLAGS) -D$(ARCH) -c $*.c
.f.o:
	$(FC) $(FFLAGS) -c $*.f
#.F.o:
#    ./Fpp $*.F $(CPPFLAGS)
#    $(FC) $(FFLAGS) -c $*.f
#    rm $*.f

cputime_$(ARCH).o: cputime.F
	./Fpp cputime.F $(CPPFLAGS)
	$(FC) $(FFLAGS) -c cputime_$(ARCH).f

