// C++ source code
// File: "/home/kassick/Work/olam/olamfs-trace-converter/include/olam_events.h"
// Created: "Qui, 02 Jun 2011 10:32:37 -0300 (kassick)"
// Updated: "Sex, 17 Jun 2011 18:11:10 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  olam_events.h
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  02-06-2011 10:32:37 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */

#ifndef _OLAM_EVENTS_H_
#define _OLAM_EVENTS_H_

#include <stdlib.h>
#include <name_types.h>
#include <macros_f.h>

#define OLAM_VERSION_BASE (10)
#define OLAM_EVT_BASE (4000)
#define OLAM_EVT_MAX  (4499)
#define EVT_BASE_NAME OLAM
#include <evt_template.h>




// Add IDs here
typedef enum _olam_evts_t {
  OLAM_INIT = OLAM_EVT_BASE,
  EVT_IO(HDF5_OPEN),
  EVT_IO(HDF5_CLOSE),
  EVT_IO(HDF5_CREATE),
  EVT_IO(HDF5_READ),
  EVT_IO(HDF5_WRITE),
  EVT_IO(THREAD),   // individual thread in/out
  EVT_IO(PARBLOCK), // whole parallel block
  EVT_IO(READ_NL),
  EVT_IO(ONAME_CHECK),
  EVT_IO(HISTORY_START),
  EVT_IO(COPY_NL),
  EVT_IO(GRIDINIT),
  EVT_IO(PARA_DECOMP),
  EVT_IO(PARA_INIT),
  EVT_IO(MODSCHED),
  EVT_IO(FILL_JTABS),
  EVT_IO(FILL_JSEA),
  EVT_IO(FILL_JLAND),
  EVT_IO(FILL_JFLUX),
  EVT_IO(ALLOC_MISC),
  EVT_IO(JNMBINIT),
  EVT_IO(O_MEM_ALLOC),
  EVT_IO(O_ALLOC_MPI),
  EVT_IO(O_ALLOC_MPI_LAND),
  EVT_IO(O_ALLOC_MPI_SEA),
  EVT_IO(INITHH),
  EVT_IO(FLDSLHI),  //WHAT!
  EVT_IO(MICINIT),
  EVT_IO(HARR_RADINIT),
  EVT_IO(FULIOU_RADINIT),
  EVT_IO(LEAF3_STARTUP),
  EVT_IO(INIT_OFFLINE_MET),
  EVT_IO(READ_OFFLINE_MET_INIT),
  EVT_IO(LEAF4_INIT_ATM),
  EVT_IO(SEA_INIT_ATM),
  EVT_IO(ISAN_DRIVER),
  EVT_IO(RAYF_INIT),
  EVT_IO(PLOTONLY),
  EVT_IO(HISTORY_WRITE),
  EVT_IO(MODEL),
  EVT_IO(O_CLSGKS),
  EVT_IO(TIMESTEP),
  EVT_IO(UPDATE_MODEL_TIME),
  EVT_IO(O_OUTPUT),
  EVT_IO(TEND0),
  EVT_IO(RADIATE),
  EVT_IO(SURFACE_TURB_FLUX),
  EVT_IO(TURB_K),
  EVT_IO(CUPARM_DRIVER),
  EVT_IO(SURFACE_CUPARM_FLUX),
  EVT_IO(THILTEND_LONG),
  EVT_IO(VELTEND_LONG),
  EVT_IO(OBS_NUDGE),
  EVT_IO(ZERO_MASSFLUX),
  EVT_IO(PROG_WRTU),
  EVT_IO(TIMEAVG_MASSFLUX),
  EVT_IO(SCALAR_TRANSPORT),
  EVT_IO(PREDTR),
  EVT_IO(THERMO),
  EVT_IO(MICRO),
  EVT_IO(SURFACE_PRECIP_FLUX),
  EVT_IO(TRSETS),
  EVT_IO(LEAF3),
  EVT_IO(SEACELLS),
  EVT_IO(INNERSTEP),




  EVT_N(TESTE),
} olam_evts_t;




// Add names to list here
static evt_name_t  olam_evt_names[] = { 
  EVT_NAME_ENTRY_N(INIT),
  EVT_NAME_ENTRY(HDF5_OPEN),
  EVT_NAME_ENTRY(HDF5_CREATE),
  EVT_NAME_ENTRY(HDF5_READ),
  EVT_NAME_ENTRY(HDF5_WRITE),
  EVT_NAME_ENTRY(THREAD),   // individual thread in/out
  EVT_NAME_ENTRY(PARBLOCK), // whole parallel block
  EVT_NAME_ENTRY(READ_NL),
  EVT_NAME_ENTRY(ONAME_CHECK),
  EVT_NAME_ENTRY(HISTORY_START),
  EVT_NAME_ENTRY(COPY_NL),
  EVT_NAME_ENTRY(GRIDINIT),
  EVT_NAME_ENTRY(PARA_DECOMP),
  EVT_NAME_ENTRY(PARA_INIT),
  EVT_NAME_ENTRY(MODSCHED),
  EVT_NAME_ENTRY(FILL_JTABS),
  EVT_NAME_ENTRY(FILL_JSEA),
  EVT_NAME_ENTRY(FILL_JLAND),
  EVT_NAME_ENTRY(FILL_JFLUX),
  EVT_NAME_ENTRY(ALLOC_MISC),
  EVT_NAME_ENTRY(JNMBINIT),
  EVT_NAME_ENTRY(O_MEM_ALLOC),
  EVT_NAME_ENTRY(O_ALLOC_MPI),
  EVT_NAME_ENTRY(O_ALLOC_MPI_LAND),
  EVT_NAME_ENTRY(O_ALLOC_MPI_SEA),
  EVT_NAME_ENTRY(INITHH),
  EVT_NAME_ENTRY(FLDSLHI),  //WHAT!
  EVT_NAME_ENTRY(MICINIT),
  EVT_NAME_ENTRY(HARR_RADINIT),
  EVT_NAME_ENTRY(FULIOU_RADINIT),
  EVT_NAME_ENTRY(LEAF3_STARTUP),
  EVT_NAME_ENTRY(INIT_OFFLINE_MET),
  EVT_NAME_ENTRY(READ_OFFLINE_MET_INIT),
  EVT_NAME_ENTRY(LEAF4_INIT_ATM),
  EVT_NAME_ENTRY(SEA_INIT_ATM),
  EVT_NAME_ENTRY(ISAN_DRIVER),
  EVT_NAME_ENTRY(RAYF_INIT),
  EVT_NAME_ENTRY(PLOTONLY),
  EVT_NAME_ENTRY(HISTORY_WRITE),
  EVT_NAME_ENTRY(MODEL),
  EVT_NAME_ENTRY(O_CLSGKS),
  EVT_NAME_ENTRY(TIMESTEP),
  EVT_NAME_ENTRY(UPDATE_MODEL_TIME),
  EVT_NAME_ENTRY(O_OUTPUT),
  EVT_NAME_ENTRY(TEND0),
  EVT_NAME_ENTRY(RADIATE),
  EVT_NAME_ENTRY(SURFACE_TURB_FLUX),
  EVT_NAME_ENTRY(TURB_K),
  EVT_NAME_ENTRY(CUPARM_DRIVER),
  EVT_NAME_ENTRY(SURFACE_CUPARM_FLUX),
  EVT_NAME_ENTRY(THILTEND_LONG),
  EVT_NAME_ENTRY(VELTEND_LONG),
  EVT_NAME_ENTRY(OBS_NUDGE),
  EVT_NAME_ENTRY(ZERO_MASSFLUX),
  EVT_NAME_ENTRY(PROG_WRTU),
  EVT_NAME_ENTRY(TIMEAVG_MASSFLUX),
  EVT_NAME_ENTRY(SCALAR_TRANSPORT),
  EVT_NAME_ENTRY(PREDTR),
  EVT_NAME_ENTRY(THERMO),
  EVT_NAME_ENTRY(MICRO),
  EVT_NAME_ENTRY(SURFACE_PRECIP_FLUX),
  EVT_NAME_ENTRY(TRSETS),
  EVT_NAME_ENTRY(LEAF3),
  EVT_NAME_ENTRY(SEACELLS),
  EVT_NAME_ENTRY(INNERSTEP),
  NULL_EVT,
};

#undef EVT_BASE_NAME

#endif
