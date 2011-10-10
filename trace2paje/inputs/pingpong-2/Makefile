# Makefile
# File: "/home/kassick/Work/olam/trace2paje/inputs/pingpong-2/Makefile"
# Created: "Qua, 05 Out 2011 14:04:42 -0300 (kassick)"
# Updated: "Seg, 10 Out 2011 17:37:19 -0300 (kassick)"
# $Id$
# Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 

BASE = pingpong


TRACE_CONVERTER=/home/kassick/Work/olam/trace2paje/src/rastro2paje
MAIN_DESC=pingpong.txt
IDS=ids.txt
EVT_FILE=pingpong-link.txt
RASTRO_GENERATE= /home/kassick/Work/librastro/bin/rastro_generate
RASTRO_INCLUDE=/home/kassick/Work/librastro/include
RASTRO_LIB=/home/kassick/Work/librastro/lib

DSC_FILE=rst_functions.dsc
EVT_DEFINES=evt_defs.h
RST_FUNCTIONS=rst_functions.c

PING_NP=8

pingpong.paje: run_pingpong
	$(TRACE_CONVERTER) -i $(MAIN_DESC) -o pingpong.paje *.rst

run_pingpong: mpi-pingpong
	mpirun -np $(PING_NP) mpi-pingpong

mpi-pingpong: mpi-pingpong.c $(RST_FUNCTIONS) $(EVT_DEFINES)
	mpicc -I$(RASTRO_INCLUDE) -L$(RASTRO_LIB) -Wl,-rpath=$(RASTRO_LIB) \
	  	-lRastro mpi-pingpong.c $(RST_FUNCTIONS)  \
		-o mpi-pingpong

$(RST_FUNCTIONS): $(DSC_FILE)
	$(RASTRO_GENERATE) -h rst_functions.h -c rst_functions.c  $(DSC_FILE)

$(DSC_FILE): $(MAIN_DESC) $(EVT_FILE) $(IDS)
	$(TRACE_CONVERTER) -i $(EVT_FILE) --rst-signatures $(DSC_FILE)

$(EVT_DEFINES): $(IDS)
	$(TRACE_CONVERTER) -i $(MAIN_DESC) --c-header $(EVT_DEFINES)

$(IDS): $(EVT_FILE)
	$(TRACE_CONVERTER) -i $(EVT_FILE) --gen-ids $(IDS)

clean:
	rm pingpong.paje rst_functions.* $(EVT_DEFINES) $(IDS) mpi-pingpong *.rst


.PHONY: clean



all: pingpong.paje
