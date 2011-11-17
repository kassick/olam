#
# Makefile for Linux machines 
#

LARGE_FILE  = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE=1 -D__USE_LARGEFILE64=1
UTILS_DIR   = ../../utilities/src
#BENCH_COMM  =-DCRBENCH_COMM=MPI_COMM_WORLD


RASTRO_DIR=/home/rkassick/Work/librastro/
#RASTRO_DIR=/usr/local
RASTRO_LIBS=-L${RASTRO_DIR}/lib
RASTRO_INCS=-I${RASTRO_DIR}/include

# get MPI_CC, MPI_INC, and MPI_LD from environment
ifdef MPI_CC 
	CC		= ${MPI_CC}	
else
	CC      = mpicc
endif
CFLAGS  = -g -lm  -lRastro -I. ${MPI_INC} ${RASTRO_INCS} ${RASTRO_LIBS} -DROMIO -DMPI_IO_TEST_21_PLUS_1 -D_GNU_SOURCE -I$(UTILS_DIR) ${LARGE_FILE} -DMPI_IO_TEST_VERSION2='"1.00.021"' -Wall
CLIBS   = ${MPI_LD}

OBJS		= burn.o
MPIIOOBJS	=	
UTILS		=	$(UTILS_DIR)/print.o $(UTILS_DIR)/utilities.o

all:    fs_test.$(MY_MPI_HOST).x

fs_test.$(MY_MPI_HOST).x: fs_test.x
	cp $< $@

fs_test.x:	$(OBJS) $(UTILS) $(MPIIOOBJS) fs_test.c
	$(CC) $(CFLAGS) $(CLIBS) $(BENCH_COMM) \
		fs_test.c \
		rastro_functs.c \
		$(MPIIOOBJS) $(UTILS) $(OBJS) $(CLIBS) \
		-o $@

%.o : %.c
	$(CC) $(CFLAGS) -c $? -o $@

$(UTILS_DIR)/%.o: $(UTILS_DIR)/%.c
	$(CC) $(CFLAGS) -c $? -o $@

clean:
	- /bin/rm -f *~
	- /bin/rm -f *.o $(UTILS_DIR)/*.o
	- /bin/rm -f fs_test.x
	- /bin/rm -f fs_test.$(MY_MPI_HOST).x
