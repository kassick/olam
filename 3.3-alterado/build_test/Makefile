.SUFFIXES:
.SUFFIXES: .F90 .f90 .o .a .c

# Include definitions

include ./paths.mk
include ./include.mk

# Compiler commands.

CC	 = icc
CFLAGS   = $(C_OPTS) $(INCLUDES) -D$(CMACH)



INCLUDES   = -I$(OUTILS) $(PAR_INCS) $(HDF5_INCS) $(NETCDF_INCS) $(RASTRO_INCS)
F_COMMAND  = $(F_COMP) -c $(F_OPTS) $(INCLUDES)
C_COMMAND  = $(C_COMP) -c $(CFLAGS)

ifndef F_DEFINE_FLAG
F_DEFINE_FLAG=-D
endif

ifdef CMACH
MACH_DEF=$(F_DEFINE_FLAG)$(CMACH)
endif

ifdef OLAM_MPI
MPI_DEF=$(F_DEFINE_FLAG)OLAM_MPI
OLAM_SUFFIX=-mpi
endif

ifdef OLAM_RASTRO
INCLUDES += $(F_DEFINE_FLAG)OLAM_RASTRO
#INCLUDES = $(INCLUDES) $(F+DEFINE_FLAG)OLAM_RASTRO
endif

# Define archive and executable names.

BASE=olam
EXE=$(BASE)-$(OLAM_VERSION)$(OLAM_SUFFIX)
ARC=$(BASE)-$(OLAM_VERSION).a

# OLAM source files

include objects.mk

# f90 rule.

.f90.a:
	@echo ""
	$(F_COMMAND) $<
	$(ARCHIVE) $@ $(<F:.f90=.o)
	rm -f $(<F:.f90=.o)

# F90 rule.

.F90.a:
	@echo ""
	$(F_COMMAND) $(MACH_DEF) $(MPI_DEF) $<
	$(ARCHIVE) $@ $(<F:.F90=.o)
	rm -f $(<F:.F90=.o)

# c rule.

.c.a:
	@echo ""
	$(C_COMMAND) $<
	$(ARCHIVE) $@ $(<F:.c=.o)
	rm -f $(<F:.c=.o)


# Define targets.

all: $(EXE)

#$(OUTILS)/rastro.c: $(OUTILS)/rastro.rsti
#	$(RASTRO_GEN) -f $(OUTILS)/rastro_f.c -c $@ -h $(OUTILS)/rastro.h     $<


MAIN_OBJ = $(notdir $(MAIN:.f90=.o))
MPIF_OBJ = $(notdir $(MPIF:.F90=.o))

$(EXE): $(ARC) $(MPIF) $(MAIN) FORCE
	@echo ""
	$(LOADER) -o $(EXE) $(MAIN_OBJ) $(MPIF_OBJ) $(LOADER_OPTS) $(ARC) \
	    $(LIBNCARG) $(PAR_LIBS) $(LIBS) $(HDF5_LIBS) $(NETCDF_LIBS) $(RASTRO_LIBS) $(OMP)
	rm -f *.o
	@echo ""
	@echo Finished building === $(EXE)
	@echo ""

$(MAIN): FORCE
	@echo ""
	$(F_COMMAND) $@

$(MPIF): FORCE
	@echo ""
	$(F_COMMAND) $(MACH_DEF) $(MPI_DEF) $@

$(ARC): $(OBJ)

FORCE:

clean:
	@echo ""
	rm -f $(ARC) $(BASE)-$(OLAM_VERSION) $(BASE)-$(OLAM_VERSION)-mpi
	rm -f *.o *.mod *.f *.stb
	@echo ""

# Include dependencies (if this does not work, try using gnu make or copy the
# contents of the file dep_rams here).

include dep_olam.mk
