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

oscila2.x: nrtype.o parametros_oscila.o nrutil_DP.o nr_DP.o \
           rkck_DP.o rkqs_DP.o odeint_DP.o oscila2.o 
	$(IFC) $(OPTION) -o oscila2.x nrtype.o parametros_oscila.o nrutil_DP.o nr_DP.o \
		                      rkck_DP.o rkqs_DP.o odeint_DP.o oscila2.o 

nrtype.o:  nrtype.f90 
	$(IFC) $(OPTION) -c   nrtype.f90 


parametros.o:  parametros.f90 
	$(IFC) $(OPTION) -c  parametros.f90 

parametros_oscila.o:  parametros_oscila.f90 
	$(IFC) $(OPTION) -c  parametros_oscila.f90 

nrutil_DP.o:  nrutil_DP.f90 
	$(IFC) $(OPTION) -c nrutil_DP.f90 

nr_DP.o:  nr_DP.f90 
	$(IFC) $(OPTION) -c nr_DP.f90 

rkck_DP.o:  rkck_DP.f90 
	$(IFC) $(OPTION) -c rkck_DP.f90 

rkqs_DP.o:  rkqs_DP.f90 
	$(IFC) $(OPTION) -c  rkqs_DP.f90 

odeint_DP.o:  odeint_DP.f90 
	$(IFC) $(OPTION) -c  odeint_DP.f90 

oscila2.o: oscila2.f90
	$(IFC) $(OPTION) -c oscila2.f90

clean:
	rm *.mod *.o *.x

clean_tmp:
	 rm *.*~ 
