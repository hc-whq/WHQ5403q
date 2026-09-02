# Makefile
#
.SUFFIXES:
.SUFFIXES: .o .f90

include ../macros

OBJS = \
	module_SedCNPvariables.o \
	module_AlloDeallocate.o \
	module_overlandSed.o \
	module_channelSed.o \
	module_WriteCNPini.o \
	module_CNPfunctions.o \
	module_CNPparams.o \
	module_SedCNP_in.o \
        module_SedCNP_out.o \
	module_Carbon.o \
	module_Nitrogen.o \
	module_Phosphorus.o \
	module_CNPupdates.o \
	module_ReadCNPini.o \
	module_CNPmain.o \
        module_SedCNPmodel_driver.o
        	

all:    $(OBJS)

.f90.o:

	@echo ""
	$(COMPILER90) $(CPPINVOKE) $(CPPFLAGS) -o $(@) $(F90FLAGS) $(LDFLAGS) $(MODFLAG) -I$(NETCDFINC) $(*).f90
	@echo ""
	ar -r ../lib/libHYDRO.a $(@)
	cp *.mod ../mod

#
# Dependencies:
#

clean:
	rm -f *.o *.mod *.stb *~

