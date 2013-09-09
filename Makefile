#OPTION=  -fpe0
SYSTEM := $(shell hostname)
#SYSTEM=debug
ifeq ($(SYSTEM),"debug")
IFC = ifort
OPTION=  -fpe0 -g 
endif

ifeq ($(SYSTEM),IA32)
IFC = ifort
endif
ifeq ($(SYSTEM),AMD)
IFC = f95i
endif
ifeq ($(SYSTEM),wepwawet)
IFC = gfortran
endif

OBJS =  nrtype.o            \
	parametros_oscila.o \
	nrutil_DP.o         \
	nr_DP.o             \
        rkck_DP.o           \
        rkqs_DP.o           \
	odeint_DP.o

oscila2.x: $(OBJS) oscila2.o 
	$(IFC) $(OPTION) -o oscila2.x $(OBJS) oscila2.o 

.SUFFIXES: $(SUFFIXES) .f90


.f90.o:
	$(IFC) $(OPTION) -c $<

clean:
	rm *.mod *.o *.x

clean_tmp:
	 rm *.*~ 
