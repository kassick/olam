// C++ source code
// File: "/home/kassick/Work/olam/olamfs-trace-converter/include/olam_events.h"
// Created: "Qui, 02 Jun 2011 10:32:37 -0300 (kassick)"
// Updated: "Qua, 08 Jun 2011 19:09:18 -0300 (kassick)"
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
  OLAM_EVT_INIT = OLAM_EVT_BASE,
  EVT_IO(HDF5_OPEN),
  EVT_IO(HDF5_CLOSE),
  EVT_IO(HDF5_CREATE),
  EVT_IO(HDF5_READ),
  EVT_IO(HDF5_WRITE),
  EVT_IO(THREAD),   // individual thread in/out
  EVT_IO(PARBLOCK), // whole parallel block
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
  NULL_EVT,
};

#undef EVT_BASE_NAME

#endif
