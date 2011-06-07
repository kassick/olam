// C++ source code
// File: "/home/kassick/Work/olam/olamfs-trace-converter/include/olam_events.h"
// Created: "Qui, 02 Jun 2011 10:32:37 -0300 (kassick)"
// Updated: "Ter, 07 Jun 2011 18:35:19 -0300 (kassick)"
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

#define OLAM_VERSION_BASE (10)
#define OLAM_EVT_BASE (4000)

// IN / OUT Declaration
#define OLAM_EVT_IO(name) \
 OLAM_EVT_ ## name ## _IN, \
 OLAM_EVT_ ## name ## _OUT

#define OLAM_EVT_N(name) \
  OLAM_EVT_ ## name 

#define OLAM_EVT_STR(id) \
  {OLAM_EVT_ ## id ## _IN, "OLAM_" #id}, \
  {OLAM_EVT_ ## id ## _OUT, "OLAM_" #id "V"}

#define OLAM_EVT_STR_N(id) \
  {OLAM_EVT_ ## id  , "OLAM_" #id}




// Add IDs here
typedef enum _olam_evts_t {
  OLAM_EVT_INIT = OLAM_EVT_BASE,
  OLAM_EVT_IO(HDF5_OPEN),
  OLAM_EVT_IO(HDF5_CLOSE),
  OLAM_EVT_IO(HDF5_CREATE),
} olam_evts_t;




typedef struct _olam_evt_names_t {
  olam_evts_t id;
  char * name;
} olam_evt_names_t;



// Add names to list here
static olam_evt_names_t  olam_evt_names[] = { 
  OLAM_EVT_STR_N(INIT),
  OLAM_EVT_STR(HDF5_OPEN),
  OLAM_EVT_STR(HDF5_CREATE),
  {0, NULL},
};


#endif
