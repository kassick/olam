# Activate appropriate parts below, comment out others:


#-----------------  LINUX Portland Group pgf77/gcc ---------------
#CMACH=PC_LINUX1
#F_COMP=pgf90
#
## If using MPI libraries:
##MPI_PATH=/usr/local/mpich
##PAR_INCS=-I$(MPI_PATH)/include
##PAR_LIBS=-L$(MPI_PATH)/lib -lmpich -lpmpich
##OLAM_MPI=yes
#
## OPTIMIZED:
#F_OPTS=-Mvect=cachesize:524288 -Munroll -Mnoframe -O2 -pc 64 \
#       -Mfree  -Mbyteswapio
#
## DEBUG:
##F_OPTS= -g  -pc 64 -Mfree  -Mbyteswapio -Mbounds
#
#C_COMP=gcc
#C_OPTS=-O3 -DUNDERSCORE -DLITTLE
#
#NCARG_DIR=/usr/local/ncarg/lib
#LIBNCARG=-L$(NCARG_DIR) -lncarg -lncarg_gks -lncarg_c \
#          -L/usr/X11/lib64 -lX11 -ldl -lpthread
#
#HDF5_LIBS=-L/usr/local/lib64 -lhdf5 -lz -lm
#HDF5_INCS=-I/usr/local/include
#
#NETCDF_LIBS=-L/usr/local/lib64 -lnetcdf
#NETCDF_INCS=-I/usr/local/include
#
#LOADER=$(F_COMP)
#LOADER_OPTS=-v -Wl,-static $(F_OPTS)
#LIBS=
#ARCHIVE=ar rs


DAMARIS=yes

ifdef DAMARIS

DAMARIS_DIR=/home/rkassick/Work/damaris-intel
DAMARIS_LIBS=-Xlinker --start-group $(DAMARIS_DIR)/lib/libdamaris.a $(DAMARIS_DIR)/lib/libdamaris-server.a -Xlinker --end-group -lboost_filesystem -lboost_system -lboost_program_options -lxerces-c -lrt -ldl
DAMARIS_INC=-I$(DAMARIS_DIR)/include $(F_DEFINE_FLAG)DAMARIS
DAMARIS_MOD=
#$(DAMARIS_DIR)/lib/damaris.mod

endif

#-----------------  LINUX Intel Fortran ifort/gcc ---------------
CMACH=PC_LINUX1
F_COMP=ifort

# If using MPI libraries:
#MPI_PATH=/usr/local/mpich
#PAR_INCS=-I/usr/include/mpich2/
PAR_INCS=-I/usr/lib/openmpi/include
#PAR_LIBS=-L/usr/lib/mpich/lib/shared/ -lmpich
#PAR_LIBS=-lmpi_f90
PAR_LIBS=-lmpichf90 -lgfortran
OLAM_MPI=yes
#OLAM_RASTRO=no

# OPTIMIZED:
F_OPTS=-fpp -xW -O3 -fno-alias -prec-div -free -traceback

# DEBUG:
#F_OPTS=-g -free -prec-div -fp -check bounds -inline-debug-info -traceback \
#        -debug extended -check uninit -ftrapuv -auto 

C_COMP=icc
C_OPTS=-O3 -DUNDERSCORE -DLITTLE

#NCARG_DIR=/usr/local/ncarg/lib
#LIBNCARG=-L$(NCARG_DIR) -lncarg -lncarg_gks -lncarg_c \
#          -L/usr/X11/lib64 -lX11 -ldl -lpthread
NCARG_DIR=/home/rkassick/Work/ncarg/lib/
LIBNCARG=-L$(NCARG_DIR) -lncarg -lncarg_gks -lncarg_c \
          -L/usr/X11/lib64 -lX11 -ldl -lpthread

         

#HDF5_LIBS=-L/usr/local/lib64 -lhdf5 -lz -lm -lgfortran
#HDF5_INCS=-I/usr/local/include

HDF5_DIR=/home/rkassick/Work/hdf5-serial-intel/
#HDF5_LIBS=-L$(HDF5_DIR)/lib -lhdf5 -lz -lm #/usr/lib/libgfortran.so.3 -- this was for hdf5 bundled with debian
HDF5_LIBS=-L$(HDF5_DIR)/lib 		  \
      -lhdf5_fortran -lhdf5		  \
      -lhdf5_hl -lhdf5hl_fortran	  \
      -lz -lm
HDF5_INCS=-I$(HDF5_DIR)/include

#NETCDF_LIBS=-L/usr/local/lib64 -lnetcdf
#NETCDF_INCS=-I/usr/local/include
NETCDF_DIR=/home/rkassick/Work/netcdf-serial-intel
NETCDF_LIBS=-L$(NETCDF_DIR)/lib -lnetcdf
NETCDF_INCS=-I$(NETCDF_DIR)/include

OMPI_MPIFC=ifort
LOADER=mpif90
LOADER_OPTS=-i-static $(F_OPTS)

# to allow ifort compiler to link with pg-compiled ncar graphics:
# LIBS=-z muldefs -L/opt/pgi/linux86-64/5.2/lib -lpgftnrtl -lpgc
ARCHIVE=ar rs

OMP=-openmp

ifdef OLAM_RASTRO
  RASTRO_LIBS=-L/home/rkassick/Work/rastro/lib -lRastro
  RASTRO_INCS=-I/home/rkassick/Work/rastro/include
  RASTRO_GEN=/home/rkassick/Work/rastro/bin/rastro_generate
endif


#----------------- IBM xlf/xlc ---------------------------------
#CMACH=AIX
#F_COMP=xlf90_r       # without MPI
#
## F_COMP=mpxlf90_r   # with MPI
## OLAM_MPI=yes       # with MPI
#
## OPTIMIZED:
#F_OPTS=-O3 -qarch=auto -qtune=auto -qcache=auto -q64 \
#       -qsuffix=cpp=F90:f=f90 -qhot -qunroll=yes -qalias=noaryovrlp,nopteovrlp
#
## DEBUG:
##F_OPTS=-O3 -g -C -qflttrap=zerodivide,invalid,ov,enable,nanq -qfullpath \
##       -qfloat=nans -qsuffix=cpp=F90:f=f90 -qstrict -qnosave -qsigtrap \
##       -qinitauto=FF -qtbtable=full -qwarn64
#
#F_DEFINE_FLAG=-WF,-D
#C_COMP=xlc_r
#C_OPTS=-O3 -qarch=auto -qtune=auto -qcache=auto -q64
#
#NCARG_DIR=/usr/local/apps/ncl-5.0.0/lib64/r4i4
#LIBNCARG=-L$(NCARG_DIR) -lncarg -lncarg_gks -lncarg_c -lX11
#
#HDF5_LIBS=-L/contrib/hdf5/lib64/r4i4 -lhdf5 -lsz -lz
#HDF5_INCS=-I/contrib/hdf5/include/hdf5-64
#
#NETCDF_LIBS=-L/usr/local/apps/netcdf-3.6.2/lib -lnetcdf
#NETCDF_INCS=-I/usr/local/apps/netcdf-3.6.2/include
#
#LOADER=$(F_COMP)
#LOADER_OPTS=$(F_OPTS)
#LIBS=
#ARCHIVE=ar rs
#
#-----------------------------------------------------------
