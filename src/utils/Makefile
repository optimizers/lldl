#
# UTILS library.
#

L_ARCH   = $(ARCH)
LIB_NAME = d-$(L_ARCH).a

CPPFLAGS  = -E -C -D$(ARCH)

AR = ar -rcv

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

