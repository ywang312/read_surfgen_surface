
OBJECTS =  evalsurf.o main.o
LIBSURF =  hddata.o combinatorial.o progdata.o libutil.o libsym.o libinternal.o \
	   io.o potlib.o getver.o
LDFLAGS =  $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
	   $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a \
	   -Wl,--start-group \
           $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
           $(MKLROOT)/lib/intel64/libmkl_core.a \
	   $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
	   -Wl,--end-group -lpthread -lm

FC:=ifort
#FC:=gfortran
#FCOPT:=-g -fdefault-integer-8 -fopenmp
FCOPT:=-O3 -qopenmp -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include

all:	lib $(OBJECTS)
	$(FC) $(FCOPT) -c evalsurf.f90
	$(FC) $(FCOPT) -o evalsurf.x evalsurf.o libsurfgen.a $(LDFLAGS)
	$(FC) $(FCOPT) -c main.f90
	$(FC) $(FCOPT) -o main.x main.o libsurfgen.a $(LDFLAGS)

lib:	$(LIBSURF)
	$(AR) -r -v libsurfgen.a $(LIBSURF)
	@echo " Cleaning up..."
	rm -f $(LIBSURF)
	rm -f ./*.mod

./%.o : ./lib/%.f90
	@echo ' Building file: $<'
	$(FC) $(FCOPT) -c -o $@ $<
	@echo ' Finished building $<'
	@echo ''

./%.o : ./lib/%.F90
	@echo ' Building file: $<'
	$(FC) $(FCOPT) -c -o $@ $<
	@echo ' Finished building $<'
	@echo ''

./%.o : ./%.f90
	@echo ' Building file: $<'
	$(FC) $(FCOPT) -c -o $@ $<
	@echo ' Finished building $<'
	@echo ''

clean :
	rm -f ./*.mod ./*.o
